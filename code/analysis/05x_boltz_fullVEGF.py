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
05x_boltz_fullVEGF.py — 全长 VEGF165 复核（关键升级任务）

目的: 之前所有 Boltz/MM-GBSA 都在 VEGF-A 1–98 截短体(受体结合核心)上做；
      HBD(99–165) 不在模型内。本脚本把「全长 VEGF165 成熟序列(165aa)」作为
      chain A，对三条已优选候选重新跑 Boltz，看全长上下文(含 HBD 存在)下
      结合预测(ipTM/界面pAE/界面pLDDT)与 98aa 截短体相比是否漂移。

候选 (chain B):
  1✱ = #1_T0.3_s132   (80aa, Flt-1+二聚面, 去A最强)
  2✱ = #2c_T0.1_s18   (100aa, KDR+Flt-1 双受体, A85E magic)
  4✱ = T0.1_s35       (80aa, KDR 单受体, 去淀粉样锁定)

输入:
  data/VEGF165_full.fasta         (165aa 全长成熟序列, 已验证 [:98]==旧靶点)
  data/binders_fullVEGF.json      (三条 binder 序列)
输出:
  results/deA_redesign_1_2/fullVEGF_boltz_3seed.csv
  results/deA_redesign_1_2/boltz_pdb_fullVEGF/<stem>/{inputs,out}/

关键补丁(继承 05k):
  --model boltz1 --no_kernels   (RTX5090 sm_120 必加)
  解析用 os.walk 避开 stem 内 [] 字符类陷阱

用法: python 05x_boltz_fullVEGF.py [--seeds 0,11,42] [--limit N]
"""
import os, sys, subprocess, glob, csv, shutil, json, argparse

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, ".."))
BOLTZ_PYTHON = r"C:\Users\Administrator\.conda\envs\boltz_cuda\python.exe"
BOLTZ_CLI = os.path.join(os.path.dirname(BOLTZ_PYTHON), "Scripts", "boltz.exe")

VEGF_FA = os.path.join(PROJ, "data", "VEGF165_full.fasta")
BINDERS_JSON = os.path.join(PROJ, "data", "binders_fullVEGF.json")
OUT_CSV = os.path.join(PROJ, "results", "deA_redesign_1_2", "fullVEGF_boltz_3seed.csv")
WORK_ROOT = os.path.join(PROJ, "results", "deA_redesign_1_2", "boltz_pdb_fullVEGF")

# 基线(98aa 截短体)结果, 用于报告中对比 (来自历史 CSV)
BASELINE_98 = {
    "1star": dict(id="#1_T0.3_s132", iptm=0.937, pae=9.88, plddt=85.58),
    "2star": dict(id="#2c_T0.1_s18", iptm=0.628, pae=8.20, plddt=79.40),
    "4star": dict(id="T0.1_s35",     iptm=0.699, pae=10.23, plddt=88.50),
}


def load_fasta_head_seq(path):
    seq = ""
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            continue
        seq += line
    return seq


def load_binders(path):
    d = json.load(open(path))
    # 固定顺序
    order = ["1star", "2star", "4star"]
    out = []
    for k in order:
        if k in d:
            out.append((k, d[k]))
    return out


def build_env():
    envdir = os.path.dirname(BOLTZ_PYTHON)
    lib_bin = os.path.join(envdir, "Library", "bin")
    scripts = os.path.join(envdir, "Scripts")
    mingw = os.path.join(envdir, "Library", "mingw-w64", "bin")
    sys32, sysroot = r"C:\Windows\System32", r"C:\Windows"
    pp = [p for p in [envdir, scripts, lib_bin, mingw, sys32, sysroot] if os.path.isdir(p)]
    return {
        "PATH": ";".join(pp), "SYSTEMROOT": sysroot, "WINDIR": sysroot,
        "SYSTEMDRIVE": "C:", "COMSPEC": os.path.join(sys32, "cmd.exe"),
        "TEMP": os.environ.get("TEMP", r"C:\Windows\Temp"),
        "TMP": os.environ.get("TMP", r"C:\Windows\Temp"),
        "USERPROFILE": os.environ.get("USERPROFILE", r"C:\Users\Administrator"),
        "HOMEDRIVE": "C:", "HOMEPATH": r"\Users\Administrator",
        "USERNAME": os.environ.get("USERNAME", "Administrator"),
        "OS": "Windows_NT", "PROCESSOR_ARCHITECTURE": "AMD64",
        "NUMBER_OF_PROCESSORS": os.environ.get("NUMBER_OF_PROCESSORS", "8"),
        "APPDATA": os.environ.get("APPDATA", ""),
        "PYTHONIOENCODING": "utf-8", "CUDA_VISIBLE_DEVICES": "0",
        "KMP_DUPLICATE_LIB_OK": "TRUE",
        "HF_ENDPOINT": "https://hf-mirror.com",
    }


def write_yaml(path, target_seq, binder_seq, stem):
    with open(path, "w", encoding="utf-8") as f:
        f.write("version: 1\n")
        f.write("sequences:\n")
        f.write("  - protein:\n")
        f.write("      id: A\n")
        f.write(f"      sequence: {target_seq}\n")
        f.write("      msa: empty\n")
        f.write("  - protein:\n")
        f.write("      id: B\n")
        f.write(f"      sequence: {binder_seq}\n")
        f.write("      msa: empty\n")


def find_confidence_json(out_dir, stem):
    base = f"confidence_{stem}_model_0.json"
    if os.path.isdir(out_dir):
        for root, _, files in os.walk(out_dir):
            for f in files:
                if f == base:
                    return os.path.join(root, f)
    hits = glob.glob(os.path.join(out_dir, "**", base), recursive=True)
    return sorted(hits)[0] if hits else None


def find_pae_npz(out_dir, stem):
    base = f"pae_{stem}_model_0.npz"
    if os.path.isdir(out_dir):
        for root, _, files in os.walk(out_dir):
            for f in files:
                if f == base:
                    return os.path.join(root, f)
    return None


def interface_pae_from_npz(npz_path, len_a):
    try:
        import numpy as np
        data = np.load(npz_path)
        pae = data[data.files[0]]
        n = pae.shape[0]
        if len_a <= 0 or len_a >= n:
            return round(float(pae.mean()), 2)
        ab = pae[:len_a, len_a:]
        ba = pae[len_a:, :len_a]
        inter = (ab.sum() + ba.sum()) / (ab.size + ba.size)
        return round(float(inter), 2)
    except Exception:
        return "NA"


def norm100(v):
    if isinstance(v, (int, float)) and v <= 1.0:
        return round(v * 100.0, 2)
    return round(v, 2) if isinstance(v, (int, float)) else v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,11,42")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default=OUT_CSV)
    args = ap.parse_args()

    seeds = [int(x) for x in args.seeds.split(",") if x != ""]
    TARGET = load_fasta_head_seq(VEGF_FA)
    LEN_A = len(TARGET)
    binders = load_binders(BINDERS_JSON)
    if args.limit:
        binders = binders[:args.limit]

    print(f"[信息] 全长 VEGF165 = {LEN_A}aa | 候选 {len(binders)} 条 | seeds={seeds}")
    # 标准全长(UniProt P15692-4): 真实1-11=APMAEGGGQNH, 12-109=旧RBD靶点(1FLT W), 110-165=HBD
    OLD_RBD = "HEVVKFMDVYQRSYCHPIETLVDIFQEYPDEIEYIFKPSCVPLMRCGGCCNDEGLECVPTEESNITMQIMRIKPHQGQHIGEMSFLQHNKCECRPKKD"
    rbd_ok = TARGET[11:109] == OLD_RBD
    n_ok = (len(TARGET) == 165) and (TARGET[:11] == "APMAEGGGQNH") and (TARGET[114] == "N")
    c_ok = TARGET[-13:] == "LNERTCRCDKPRR"
    print(f"[信息] 标准全长校验: 长度165={len(TARGET)==165} | RBD段(12-109)==旧靶点:{rbd_ok} | 115位=N:{n_ok} | C端完整:{c_ok}")
    assert rbd_ok and n_ok and c_ok, "全长序列校验失败，中止!"

    env = build_env()
    os.makedirs(WORK_ROOT, exist_ok=True)
    rows = []

    for si, (tag, seq) in enumerate(binders, 1):
        seed_metrics = []
        for seed in seeds:
            stem = f"{tag}_fullVEGF_s{seed}"
            in_dir = os.path.join(WORK_ROOT, stem, "inputs")
            out_dir = os.path.join(WORK_ROOT, stem, "out")
            os.makedirs(in_dir, exist_ok=True)
            yaml_path = os.path.join(in_dir, f"{stem}.yaml")
            write_yaml(yaml_path, TARGET, seq, stem)

            cmd = [BOLTZ_CLI, "predict", in_dir, "--out_dir", out_dir,
                   "--output_format", "pdb", "--recycling_steps", "3",
                   "--write_full_pae", "--model", "boltz1", "--no_kernels",
                   "--override", "--seed", str(seed)]
            print(f"[{si}/{len(binders)}] {tag} seed={seed}  boltz predict (A={LEN_A}aa B={len(seq)}aa)", flush=True)
            r = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=1800)
            if r.returncode != 0:
                print(f"  FAIL seed={seed}: {r.stderr[-600:]}", flush=True)
                seed_metrics.append((seed, None, None, None, None))
                continue
            cj = find_confidence_json(out_dir, stem)
            iptm = ptm = iplddt = None
            if cj:
                d = json.load(open(cj))
                iptm = d.get("iptm"); ptm = d.get("ptm")
                iplddt = d.get("complex_iplddt") if d.get("complex_iplddt") is not None else d.get("complex_plddt")
            ipae = "NA"
            npz = find_pae_npz(out_dir, stem)
            if npz:
                ipae = interface_pae_from_npz(npz, LEN_A)
            pdbs = glob.glob(os.path.join(out_dir, "**", f"{stem}_model_0.pdb"), recursive=True)
            if pdbs:
                shutil.copy2(pdbs[0], os.path.join(WORK_ROOT, os.path.basename(pdbs[0])))
            iplddt = norm100(iplddt) if iplddt is not None else None
            seed_metrics.append((seed, iptm, ptm, ipae, iplddt))
            print(f"  -> seed={seed} ipTM={iptm} ptm={ptm} 界面pAE={ipae} 界面pLDDT={iplddt}", flush=True)

        def avg(vals):
            vs = [v for v in vals if isinstance(v, (int, float))]
            return round(sum(vs) / len(vs), 3) if vs else None
        iptms = [m[1] for m in seed_metrics]
        paes = [m[3] for m in seed_metrics if isinstance(m[3], (int, float))]
        pls = [m[4] for m in seed_metrics if isinstance(m[4], (int, float))]
        m_iptm = avg(iptms); m_pae = avg(paes); m_pl = avg(pls)
        base = BASELINE_98.get(tag, {})
        rows.append({
            "id": tag, "len": len(seq), "lenA": LEN_A,
            "iptm_mean": m_iptm if m_iptm is not None else "NA",
            "iptm_best": max([v for v in iptms if isinstance(v, (int, float))], default="NA"),
            "ptm_mean": avg([m[2] for m in seed_metrics]),
            "interface_pae_mean": m_pae if m_pae is not None else "NA",
            "interface_plddt_mean": m_pl if m_pl is not None else "NA",
            "n_seed": len(seeds),
            "base98_iptm": base.get("iptm", "NA"),
            "delta_iptm": round(m_iptm - base["iptm"], 3) if (m_iptm is not None and "iptm" in base) else "NA",
        })

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["id", "len", "lenA", "iptm_mean", "iptm_best",
                              "ptm_mean", "interface_pae_mean", "interface_plddt_mean",
                              "n_seed", "base98_iptm", "delta_iptm"])
        w.writeheader()
        w.writerows(rows)
    print(f"\n[成功] 全长 VEGF165 Boltz -> {args.out}")
    for r in rows:
        print(f"  {r['id']:8s} full_ipTM={r['iptm_mean']} (98aa={r['base98_iptm']}, Δ={r['delta_iptm']}) "
              f"界面pAE={r['interface_pae_mean']} 界面pLDDT={r['interface_plddt_mean']}")


if __name__ == "__main__":
    main()
