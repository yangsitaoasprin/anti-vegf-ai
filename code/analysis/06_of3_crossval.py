# =============================================================================
# Author : Yang Sitao (杨四涛)
# Affiliation : College of Pharmacy, Dali University;
#               The Third Affiliated Hospital of Dali University
# Contact : yangsitaoasprin@swpu.edu.cn
# License : MIT (see ../../LICENSE)
# =============================================================================
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
06_of3_crossval.py — OpenFold3 桥接：与 Boltz-1 独立交叉验证 + VEGF 二聚体建模

设计目标
--------
我们整套从头设计流程一直只用 Boltz-1 单方法，且 RTX5090 seed 不稳定反复产出穿模伪结构。
OpenFold3（AF3 开源复现，部署在 WSL Ubuntu-24.04 / conda env openfold）作为**完全独立**的
预测器，可提供：
  - chain_pair_iptm：逐链对界面分，单独抽出 VEGF<->binder 界面
  - has_clash：内置穿模标志（0/1），与我们 MM-GBSA vdW 门控互补
  - iptm / ptm / avg_plddt：全局指标
这构成反驳“伪结构/过拟合”的独立佐证；并可用多链建模直接验证 VEGF 二聚阻断假设。

调用链
------
Windows 侧生成 AF3 风格 query JSON -> wsl -d Ubuntu-24.04 调 run_openfold predict
-> 解析 <out>/<query>/seed_*/<query>_seed_*_confidences_aggregated.json -> 汇总 CSV。

注意
----
- of3-p2-155k.pt 是 155k 步训练检查点（非最终 AF3 等效权重），OF3 结果作**佐证**，
  最终结合能裁决仍走 MM-GBSA vdW 门控。MSA 默认单序列（--use-msa-server=false）。
- OF3 模板机制（比 AF3 简单）：每条 chain 加 `template_cif_paths` + `template_cif_chain_ids`
  两个字段，OF3 自动按序列对齐；无需手填残基索引。本脚本 --use-template 即对 VEGF 链注入
  实验模板（默认 data/1FLT.cif，链 A = VEGF RBD），用于抬 VEGF 单体建模质量、补偿 ipTM 偏低。
- --mode ablation：对 1✱ 做「单序列 / 开模板 / 开 MSA」三条件对照，量化模板与 MSA 对 ipTM 的提升。
"""
import os, json, glob, argparse, subprocess, sys, csv

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# 标准成熟 VEGF-A165 (UniProt P15692-4 isoform; 三源交叉验证: Abcam/Cusabio/UniProt)
# 真实1-11=APMAEGGGQNH, 12-109=旧RBD靶点(1FLT W), 110-165=HBD(115位=单N, 非206插入段)
VEGF165 = ("APMAEGGGQNHHEVVKFMDVYQRSYCHPIETLVDIFQEYPDEIEYIFKPSCVPLMRCGGCCNDEGLECVPTEE"
           "SNITMQIMRIKPHQGQHIGEMSFLQHNKCECRPKKDRARQENPCGPCSERRKHLFVQDPQTCKCSCKNTDSRCK"
           "ARQLELNERTCRCDKPRR")  # 165 aa 标准全长成熟序列 (2026-08-18 修正)
BINDERS = {
    "1star": "GPGIALCRAAVARLRAAADAIAAAGSVDAAFAAADALPSDFAPLAQALRANAAAGTPAETAAALRDAALALEKECQAMEA",
    "2star": "LEALKEQLLKYAKTAEKCAKEAQAKAEEEQKKSEYYKKIAEDEKLAEEFAKEIGGISVEEAKEIAKDLAEAYKLAAEHYKRTSEECQARAEALKAQAAAL",
    "4star": "MSTEQAKAHTEAIEKALAIAKTEGAAAGAAYAADFMATLSLEDAASFAAMAVEAARAEGLAPEDVALFAAAGEAAVAAKR",
}
CKPT = "/root/.openfold3/of3-p2-155k.pt"
WSL_DISTRO = "Ubuntu-24.04"
OUT_ROOT = os.path.join(PROJ, "results", "deA_redesign_1_2", "of3_out")
DEFAULT_TMPL_CIF = os.path.join(PROJ, "data", "1FLT.cif")
DEFAULT_TMPL_CHAIN = "A"  # 1FLT 中 VEGF 为二聚体，链 A/B 均为 VEGF（链 C/D 为 Flt-1）


def to_wsl_path(p: str) -> str:
    p = p.replace("\\", "/")
    if len(p) >= 2 and p[1] == ":":
        return "/mnt/" + p[0].lower() + p[2:]
    return p


def chain_obj(cids, mol, seq, tmpl_paths=None, tmpl_chids=None):
    d = {"molecule_type": mol, "chain_ids": list(cids), "sequence": seq}
    if tmpl_paths:
        d["template_cif_paths"] = tmpl_paths
        d["template_cif_chain_ids"] = tmpl_chids
    return d


def build_query(name: str, chains, seeds):
    return {"seeds": list(seeds), "queries": {name: {"chains": chains}}}


def crossval_specs():
    """top-3 binder + 全长 VEGF165（双链：A=VEGF, B=binder）"""
    out = []
    for bid, seq in BINDERS.items():
        name = f"{bid}_of3_cv"
        chains = [chain_obj(["A"], "protein", VEGF165), chain_obj(["B"], "protein", seq)]
        out.append((name, chains))
    return out


def dimer_specs():
    """2×VEGF165 + binder（三链：A=VEGF, B=VEGF, C=binder），测二聚阻断"""
    out = []
    for bid, seq in BINDERS.items():
        name = f"{bid}_of3_dimer"
        chains = [chain_obj(["A"], "protein", VEGF165), chain_obj(["B"], "protein", VEGF165),
                  chain_obj(["C"], "protein", seq)]
        out.append((name, chains))
    return out


def ablation_specs(tmpl_cif_wsl, tmpl_chain):
    """1✱ 三条件对照：单序列 / 开模板 / 开 MSA。返回 (cond, chains, use_msa, use_tmpl)"""
    out = []
    base = [chain_obj(["A"], "protein", VEGF165), chain_obj(["B"], "protein", BINDERS["1star"])]
    out.append(("ss", base, False, False))  # single-sequence baseline
    if tmpl_cif_wsl:
        ch = [chain_obj(["A"], "protein", VEGF165), chain_obj(["B"], "protein", BINDERS["1star"])]
        ch[0]["template_cif_paths"] = [tmpl_cif_wsl]
        ch[0]["template_cif_chain_ids"] = [tmpl_chain]
        out.append(("tmpl", ch, False, True))
    out.append(("msa", base, True, False))  # ColabFold MSA server（需联网）
    return out


def run_of3(query_json_win: str, out_dir_win: str, n_seeds: int, name: str,
            use_msa_server: bool = False, use_templates: bool = False):
    q_wsl = to_wsl_path(query_json_win)
    o_wsl = to_wsl_path(out_dir_win)
    msa = "true" if use_msa_server else "false"
    tmpl = "true" if use_templates else "false"
    inner = (f'source /root/miniconda3/etc/profile.d/conda.sh && '
             f'conda activate openfold && '
             f'run_openfold predict --query-json {q_wsl} '
             f'--inference-ckpt-path {CKPT} --use-msa-server={msa} '
             f'--use-templates={tmpl} --num-diffusion-samples 1 '
             f'--num-model-seeds {n_seeds} --output-dir {o_wsl}')
    cmd = ["wsl", "-d", WSL_DISTRO, "--", "bash", "-lc", inner]
    log = os.path.join(out_dir_win, f"{name}.of3.log")
    os.makedirs(out_dir_win, exist_ok=True)
    print(f"[run_of3] {name} msa={msa} tmpl={tmpl} n_seeds={n_seeds} log={log}", flush=True)
    try:
        with open(log, "wb") as lf:
            r = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, timeout=900)
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"run_openfold timeout (900s) for {name}")
    try:
        tail = open(log, "r", encoding="utf-8", errors="replace").read()[-2000:]
    except Exception:
        tail = ""
    if r.returncode != 0:
        sys.stderr.write(tail + "\n")
        raise RuntimeError(f"run_openfold failed rc={r.returncode} for {name}")
    return r


def parse_aggregated(out_dir_win: str, query_name: str):
    pat = os.path.join(out_dir_win, query_name, "seed_*",
                       f"{query_name}_seed_*_confidences_aggregated.json")
    files = sorted(glob.glob(pat))
    rows = []
    for f in files:
        d = json.load(open(f))
        seed_dir = os.path.basename(os.path.dirname(f))
        rows.append({
            "query": query_name,
            "seed_dir": seed_dir,
            "iptm": d.get("iptm"),
            "ptm": d.get("ptm"),
            "avg_plddt": d.get("avg_plddt"),
            "has_clash": d.get("has_clash"),
            "disorder": d.get("disorder"),
            "chain_ptm": d.get("chain_ptm"),
            "chain_pair_iptm": d.get("chain_pair_iptm"),
        })
    return rows


def _inject_template(chains, tmpl_wsl, tmpl_chain, is_dimer):
    """给 VEGF 链注入实验模板（cv: 仅 A；dimer: A、B 两条 VEGF 均注入）"""
    out = []
    for i, c in enumerate(chains):
        if c["molecule_type"] == "protein" and c["sequence"] == VEGF165:
            if is_dimer and i >= 2:
                out.append(c)  # 第 3 条是 binder，不注模板
                continue
            c = dict(c)
            c["template_cif_paths"] = [tmpl_wsl]
            c["template_cif_chain_ids"] = [tmpl_chain]
        out.append(c)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["test", "crossval", "dimer", "all", "ablation"], default="all")
    ap.add_argument("--seeds-crossval", type=int, default=3)
    ap.add_argument("--seeds-dimer", type=int, default=1)
    ap.add_argument("--out-root", default=OUT_ROOT)
    ap.add_argument("--skip-existing", action="store_true")
    ap.add_argument("--use-msa", action="store_true",
                    help="ablation 模式包含 MSA 条件（--use-msa-server=true，需联网 ColabFold）")
    ap.add_argument("--use-template", action="store_true",
                    help="常规模式也给 VEGF 链加 1FLT 模板（修复 OF3 对 VEGF 单体建模弱导致的低 ipTM）；query 名加 _tmpl 后缀避免与 ss 版冲突")
    ap.add_argument("--template-cif", default=DEFAULT_TMPL_CIF,
                    help="VEGF 实验模板 mmCIF（Windows 路径，自动转 WSL）")
    ap.add_argument("--template-chain-id", default=DEFAULT_TMPL_CHAIN,
                    help="模板 mmCIF 中 VEGF 所在链 ID（1FLT 为 A）")
    args = ap.parse_args()

    os.makedirs(args.out_root, exist_ok=True)
    tmpl_wsl = to_wsl_path(args.template_cif) if (args.use_template and os.path.exists(args.template_cif)) else None
    if args.use_template and not tmpl_wsl:
        print(f"[WARN] 模板文件不存在 {args.template_cif}，--use-template 失效", flush=True)
    all_rows = []
    specs = []

    if args.mode in ("test", "crossval", "all"):
        for name, chains in crossval_specs():
            if tmpl_wsl:
                name2 = name + "_tmpl"
                chains2 = _inject_template(chains, tmpl_wsl, args.template_chain_id, is_dimer=False)
                specs.append(("cv_tmpl", name2, chains2, args.seeds_crossval, False, True))
            else:
                specs.append(("cv", name, chains, args.seeds_crossval, False, False))
    if args.mode in ("test", "dimer", "all"):
        for name, chains in dimer_specs():
            if tmpl_wsl:
                name2 = name + "_tmpl"
                chains2 = _inject_template(chains, tmpl_wsl, args.template_chain_id, is_dimer=True)
                specs.append(("dimer_tmpl", name2, chains2, args.seeds_dimer, False, True))
            else:
                specs.append(("dimer", name, chains, args.seeds_dimer, False, False))
    if args.mode == "ablation":
        tmpl = to_wsl_path(args.template_cif) if os.path.exists(args.template_cif) else None
        if not tmpl:
            print(f"[WARN] 模板文件不存在 {args.template_cif}，跳过 tmpl 条件", flush=True)
        for cond, chains, use_msa, use_tmpl in ablation_specs(tmpl, args.template_chain_id):
            if cond == "msa" and not args.use_msa:
                continue
            specs.append(("ablation", f"1star_{cond}", chains, args.seeds_crossval, use_msa, use_tmpl))
    if args.mode == "test":
        specs = specs[:1]

    for kind, name, chains, n_seeds, use_msa, use_tmpl in specs:
        qj = os.path.join(args.out_root, f"{name}.query.json")
        od = args.out_root
        json.dump(build_query(name, chains, list(range(n_seeds))), open(qj, "w"), indent=2)
        if args.skip_existing and glob.glob(os.path.join(od, name, "seed_*", "*aggregated.json")):
            print(f"[skip] {name} 已存在", flush=True)
        else:
            try:
                run_of3(qj, od, n_seeds, name, use_msa, use_tmpl)
            except Exception as e:
                print(f"[WARN] {name} 运行失败: {e}", flush=True)
        rows = parse_aggregated(args.out_root, name)
        for r in rows:
            r["kind"] = kind
            if kind == "ablation":
                r["condition"] = name.split("_")[-1]
        all_rows += rows
        print(f"[done] {name}: {len(rows)} seed(s) -> "
              f"iptm={[round(x['iptm'],3) if x['iptm'] is not None else None for x in rows]} "
              f"has_clash={[x['has_clash'] for x in rows]}", flush=True)

    if args.mode == "ablation":
        csv_path = os.path.join(args.out_root, "of3_ablation.csv")
        cols = ["kind", "condition", "query", "seed_dir", "iptm", "ptm", "avg_plddt",
                "has_clash", "disorder", "chain_ptm", "chain_pair_iptm"]
    else:
        csv_path = os.path.join(args.out_root, "of3_crossval.csv")
        cols = ["kind", "query", "seed_dir", "iptm", "ptm", "avg_plddt",
                "has_clash", "disorder", "chain_ptm", "chain_pair_iptm"]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in all_rows:
            w.writerow(r)
    print(f"[csv] wrote {csv_path} ({len(all_rows)} rows)", flush=True)


if __name__ == "__main__":
    main()
