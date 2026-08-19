# =============================================================================
# Author : Yang Sitao (杨四涛)
# Affiliation : School of Pharmaceutical Sciences, Dali University;
#               The Third Affiliated Hospital of Dali University
# Contact : yangsitaoasprin@swpu.edu.cn
# License : MIT (see ../../LICENSE)
# =============================================================================
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
04_boltz_rerank.py
============================================================
精度层精排（AF3 替代）：用开源 Boltz-1（AlphaFold3 级模型）对 ESMFold
入围的 binder 序列 + VEGF-A 靶点做复合物结构预测，解析界面置信度指标
（ipTM / 界面 pAE / 界面 pLDDT / ptm），按阈值过滤。

为什么用 Boltz-1 替代 AF3
------------------------
AF3 权重需 DeepMind 学术许可、本机无法直接运行；Boltz-1 是开源、可本地
GPU 运行的 AF3 级复合物预测模型，产出**同口径**的 ipTM/pAE/pLDDT，
输出 CSV 与 04_af3_rerank.py 完全兼容 → 下游 05 排序无需改动。

Windows / RTX 5090 (sm_120) 关键补丁
------------------------------------
Boltz 默认 `use_kernels=True`，会调用编译型 CUDA 算子
cuequivariance_ops_torch —— 该包**没有 Windows / py3.11 的预编译 wheel**，
安装即失败。本脚本通过以下两点规避：
  1. `--model boltz1`：使用 Boltz-1 架构（无额外许可依赖）。
  2. `--no_kernels`：强制 `use_kernels=False`，改用纯 PyTorch einsum 实现
     的三角乘/三角注意力（功能完全一致，仅稍慢；RTX 5090 上小复合物足够快）。
（注：Boltz 的 setup() 仅在 GPU 算力 < 8.0 时自动关 kernel；RTX 5090 是
 sm_120=major12 ≥ 8.0，不会自动关，必须显式 --no_kernels。）

数据流
------
  02b RFdiffusion 复合物 → 02c ProteinMPNN 序列 → 03 ESMFold 粗筛(pLDDT)
    → [本脚本] Boltz-1 精排(ipTM/界面pAE) → 05 排序

单序列模式（默认）
----------------
de novo binder 无同源序列，靶点 MSA 亦非必需；YAML 显式写 `msa: empty`
（Boltz 以序列本身作为唯一 MSA 行），全本地、无需联网。加 --use-msa-server
可对靶点走在线 MMseqs2 MSA（需联网）。

用法
----
  python 04_boltz_rerank.py \
      --seqs ../designs/mpnn \
      --prescreen-csv ../screen/esmfold/prescreen_results.csv \
      --target-pdb ../data/1FLT.pdb --target-chain W \
      --out ../results/boltz_rerank.csv \
      --min-iptm 0.6 --max-ipae 10 --min-iplddt 70

  # 数据流验证（不跑模型）
  python 04_boltz_rerank.py --seqs ../designs/mpnn --demo

依赖: conda env `boltz_cuda`（由 esm 克隆而来，含 CUDA torch 2.12+cu128 +
pip install boltz），首次运行自动下载权重（已 patch 走 HF 镜像 hf-mirror.com）。
"""
import argparse
import csv
import glob
import json
import os
import site
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.normpath(os.path.join(HERE, ".."))
BOLTZ_PYTHON = os.environ.get("BOLTZ_PYTHON",
                              r"C:\Users\Administrator\.conda\envs\boltz_cuda\python.exe")

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U",
}


# ---------------------------------------------------------------------------
# 精简子进程环境（绕过 Windows 环境变量超长）
# ---------------------------------------------------------------------------
def build_env(extra=None):
    envdir = os.path.dirname(BOLTZ_PYTHON)
    lib_bin = os.path.join(envdir, "Library", "bin")
    scripts = os.path.join(envdir, "Scripts")
    mingw = os.path.join(envdir, "Library", "mingw-w64", "bin")
    sys32 = r"C:\Windows\System32"
    sysroot = r"C:\Windows"
    path_parts = [p for p in [envdir, scripts, lib_bin, mingw, sys32, sysroot]
                  if os.path.isdir(p)]
    env = {
        "PATH": ";".join(path_parts),
        "SYSTEMROOT": sysroot, "WINDIR": sysroot, "SYSTEMDRIVE": "C:",
        "COMSPEC": os.path.join(sys32, "cmd.exe"),
        "TEMP": os.environ.get("TEMP", r"C:\Windows\Temp"),
        "TMP": os.environ.get("TMP", r"C:\Windows\Temp"),
        "USERPROFILE": os.environ.get("USERPROFILE", r"C:\Users\Administrator"),
        "HOMEDRIVE": "C:", "HOMEPATH": r"\Users\Administrator",
        "USERNAME": os.environ.get("USERNAME", "Administrator"),
        "OS": "Windows_NT", "PROCESSOR_ARCHITECTURE": "AMD64",
        "NUMBER_OF_PROCESSORS": os.environ.get("NUMBER_OF_PROCESSORS", "8"),
        "PYTHONPATH": "", "PYTHONIOENCODING": "utf-8",
        "CUDA_VISIBLE_DEVICES": os.environ.get("CUDA_VISIBLE_DEVICES", "0"),
        "KMP_DUPLICATE_LIB_OK": "TRUE",
        # Boltz 权重缓存目录 + HF 镜像（国内网络）
        "HF_ENDPOINT": os.environ.get("HF_ENDPOINT", "https://hf-mirror.com"),
    }
    py_paths = []
    try:
        py_paths += list(site.getsitepackages())
    except Exception:
        pass
    try:
        us = site.getusersitepackages()
        if us:
            py_paths.append(us)
    except Exception:
        pass
    py_paths = [p for p in py_paths if os.path.isdir(p)]
    if py_paths:
        env["PYTHONPATH"] = ";".join(py_paths)
    if extra:
        env.update(extra)
    return env


# ---------------------------------------------------------------------------
# 序列读取
# ---------------------------------------------------------------------------
def seq_from_pdb(pdb_path, chain_id):
    """从 PDB 指定链提取单字母序列（按残基编号排序，仅标准氨基酸）。"""
    seen = {}
    for line in open(pdb_path):
        if line.startswith("ATOM") and len(line) > 21 and line[21] == chain_id:
            resn = line[17:20].strip()
            try:
                resi = int(line[22:26])
            except ValueError:
                continue
            if resi not in seen and resn in AA3TO1:
                seen[resi] = AA3TO1[resn]
    return "".join(seen[k] for k in sorted(seen))


def read_fastas(seq_dir):
    records = []
    for path in glob.glob(os.path.join(seq_dir, "**", "*.fa*"), recursive=True):
        header, seq = None, []
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    if header is not None:
                        records.append((header.split()[0], "".join(seq)))
                    header, seq = line[1:], []
                elif line:
                    seq.append(line)
        if header is not None:
            records.append((header.split()[0], "".join(seq)))
    return records


def load_survivors(prescreen_csv):
    """读 03 粗筛 CSV，返回 pass==1 的 id 集合；无文件则返回 None(表示全用)。"""
    if not prescreen_csv or not os.path.exists(prescreen_csv):
        return None
    keep = set()
    with open(prescreen_csv, encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if str(row.get("pass", "")).strip() in ("1", "True", "true"):
                keep.add(row.get("id", "").strip())
    return keep


# ---------------------------------------------------------------------------
# Boltz 输入 YAML
# ---------------------------------------------------------------------------
def write_boltz_yaml(path, target_seq, binder_seq, use_msa_server):
    # de novo binder 无同源序列 → 单序列模式：显式写 msa: empty，
    # 告知 Boltz 以序列本身作为唯一 MSA 行（全本地、无需 --use_msa_server、无需联网）。
    # 注意：省略 msa 键会让 Boltz 默认去自动生成 MSA 而报错。
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("version: 1\n")
        fh.write("sequences:\n")
        msa_line = "" if use_msa_server else "      msa: empty\n"
        fh.write("  - protein:\n")
        fh.write("      id: A\n")
        fh.write(f"      sequence: {target_seq}\n")
        fh.write(msa_line)
        fh.write("  - protein:\n")
        fh.write("      id: B\n")
        fh.write(f"      sequence: {binder_seq}\n")
        fh.write(msa_line)


# ---------------------------------------------------------------------------
# 解析 Boltz 输出
# ---------------------------------------------------------------------------
def find_confidence_json(out_dir, stem):
    pats = [
        os.path.join(out_dir, f"boltz_results_{stem}", "predictions", stem,
                     f"confidence_{stem}_model_0.json"),
        os.path.join(out_dir, "**", f"confidence_{stem}_model_0.json"),
        os.path.join(out_dir, "**", f"confidence_{stem}*.json"),
    ]
    for p in pats:
        hits = glob.glob(p, recursive=True)
        if hits:
            return sorted(hits)[0]
    return None


def find_pae_npz(out_dir, stem):
    hits = glob.glob(os.path.join(out_dir, "**", f"pae_{stem}_model_0.npz"),
                     recursive=True)
    return sorted(hits)[0] if hits else None


def interface_pae_from_npz(npz_path, len_a):
    """计算跨链(A-B)界面 pAE 均值。需 numpy。"""
    try:
        import numpy as np
        data = np.load(npz_path)
        pae = data[data.files[0]]  # (N,N)
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
    """pLDDT 若为 0-1 则转 0-100。"""
    if isinstance(v, (int, float)) and v <= 1.0:
        return round(v * 100.0, 2)
    return round(v, 2) if isinstance(v, (int, float)) else v


def check_boltz_weights():
    """权重完整性预检：Boltz-1 置信度权重约 3.60GB，hf-mirror 下载易截断。
    若缓存文件明显偏小，提前报错，避免跑到模型加载时才崩在 'checkpoint corrupted'。"""
    cache = os.environ.get("BOLTZ_CACHE") or os.path.expanduser("~/.boltz")
    ckpt = os.path.join(cache, "boltz1_conf.ckpt")
    if not os.path.exists(ckpt):
        print(f"[错误] 未找到 Boltz-1 权重: {ckpt}\n"
              f"        请先用 boltz predict 触发下载，或手动从 "
              f"https://hf-mirror.com/boltz-community/boltz-1/resolve/main/boltz1_conf.ckpt 下载。")
        sys.exit(3)
    size = os.path.getsize(ckpt)
    # 服务器真实大小 ~3.60GB；小于 3.0GB 视为截断
    if size < 3_000_000_000:
        print(f"[错误] Boltz-1 权重疑似下载截断: {ckpt}\n"
              f"        本地大小 {size/1e9:.2f}GB，预期 ~3.60GB。\n"
              f"        请删除该文件后重新下载（例如用 curl -L -C - 断点续传）。")
        sys.exit(3)
    print(f"[信息] Boltz-1 权重检查通过 ({size/1e9:.2f}GB)")


def main():
    ap = argparse.ArgumentParser(description="Boltz-1 精排（AF3 替代，--no_kernels 绕过 cuequivariance）")
    ap.add_argument("--seqs", required=True, help="binder 序列 FASTA 目录（02c 产出）")
    ap.add_argument("--prescreen-csv", default=None,
                    help="03 ESMFold 粗筛 CSV，仅精排 pass==1 的入围者（可选）")
    ap.add_argument("--target-pdb", default=os.path.join(PROJECT, "data", "1FLT.pdb"))
    ap.add_argument("--target-chain", default="W", help="靶点链（1FLT 中 VEGF-A = W）")
    ap.add_argument("--target-seq", default=None, help="直接给靶点序列（优先于 PDB 提取）")
    ap.add_argument("--out", default=os.path.join(PROJECT, "results", "boltz_rerank.csv"))
    ap.add_argument("--min-iptm", type=float, default=0.6)
    ap.add_argument("--max-ipae", type=float, default=10.0)
    ap.add_argument("--min-iplddt", type=float, default=70.0)
    ap.add_argument("--use-msa-server", action="store_true",
                    help="靶点走在线 MMseqs2 MSA（需联网）；默认单序列全本地")
    ap.add_argument("--recycling", type=int, default=3)
    ap.add_argument("--limit", type=int, default=0, help="最多精排前 N 条（0=全部）")
    ap.add_argument("--demo", action="store_true")
    args = ap.parse_args()

    # 靶点序列
    if args.target_seq:
        target_seq = args.target_seq.strip()
    else:
        target_seq = seq_from_pdb(args.target_pdb, args.target_chain)
    print(f"[信息] 靶点 VEGF-A 序列长度 = {len(target_seq)}")

    # binder 序列
    records = read_fastas(args.seqs)
    survivors = load_survivors(args.prescreen_csv)
    if survivors is not None:
        records = [(i, s) for (i, s) in records if i in survivors]
        print(f"[信息] 按粗筛入围过滤后 binder = {len(records)} 条")
    else:
        print(f"[信息] 未提供粗筛 CSV，全部 binder = {len(records)} 条")
    if args.limit > 0:
        records = records[:args.limit]

    if args.demo:
        print("[DEMO] 数据流验证（不跑 Boltz）:")
        for i, (sid, seq) in enumerate(records[:5]):
            print(f"  {sid}: 复合物 A({len(target_seq)}) + B({len(seq)}) → ipTM/pAE/pLDDT")
        print(f"[DEMO] 将输出 → {args.out}（schema 同 04_af3_rerank.py）")
        return

    if not records:
        print("[错误] 无可精排的 binder 序列"); sys.exit(1)
    if not target_seq:
        print("[错误] 靶点序列为空，检查 --target-pdb/--target-chain"); sys.exit(1)

    # 权重完整性预检（避免下载截断导致 'checkpoint corrupted'）
    check_boltz_weights()

    # 写所有 YAML 到临时目录，一次性交 Boltz 批量预测
    work = tempfile.mkdtemp(prefix="boltz_")
    in_dir = os.path.join(work, "inputs")
    out_dir = os.path.join(work, "out")
    os.makedirs(in_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    id_map = {}
    for sid, seq in records:
        if not seq or any(c not in "ACDEFGHIKLMNPQRSTVWY" for c in seq):
            continue
        stem = sid.replace(".", "_")
        write_boltz_yaml(os.path.join(in_dir, f"{stem}.yaml"),
                         target_seq, seq, args.use_msa_server)
        id_map[stem] = (sid, len(seq))

    print(f"[信息] 写出 {len(id_map)} 个 Boltz 输入，开始预测（GPU）...")
    env = build_env()
    boltz_cli = os.path.join(os.path.dirname(BOLTZ_PYTHON), "Scripts", "boltz.exe")
    cmd = [boltz_cli, "predict", in_dir,
           "--out_dir", out_dir,
           "--output_format", "pdb",
           "--recycling_steps", str(args.recycling),
           "--write_full_pae",
           "--model", "boltz1",     # AF3 替代：Boltz-1 架构
           "--no_kernels",          # 强制纯 PyTorch，绕过 cuequivariance_ops_torch（无 Windows wheel）
           "--override"]
    if args.use_msa_server:
        cmd.append("--use_msa_server")
    print("  $ " + " ".join(cmd), flush=True)

    # 把 Boltz 的完整 stdout/stderr 落盘，便于排查（之前被 capture_output 吞掉）
    log_path = args.out + ".boltz.log"
    os.makedirs(os.path.dirname(log_path) or ".", exist_ok=True)
    with open(log_path, "w", encoding="utf-8") as lf:
        lf.write("$ " + " ".join(cmd) + "\n\n")
        r = subprocess.run(cmd, env=env, stdout=lf, stderr=subprocess.STDOUT, text=True)
    if r.returncode != 0:
        print(f"[错误] Boltz 预测失败（returncode={r.returncode}）。完整日志见: {log_path}")
        # 把关键报错尾部打到控制台
        try:
            with open(log_path, encoding="utf-8") as lf:
                tail = lf.read()[-2500:]
            print("---- Boltz 日志尾部 ----\n" + tail)
        except Exception:
            pass
        sys.exit(2)
    print(f"[信息] Boltz 预测完成，日志: {log_path}")

    # 解析每个预测
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    rows = []
    for stem, (sid, blen) in id_map.items():
        cj = find_confidence_json(out_dir, stem)
        if not cj:
            rows.append({"id": sid, "iptm": "NA", "ptm": "NA", "ranking": "NA",
                         "interface_pae": "NA", "interface_plddt": "NA",
                         "buried_area_A2": "NA", "pass": 0})
            continue
        with open(cj, encoding="utf-8") as fh:
            c = json.load(fh)
        iptm = c.get("iptm", c.get("protein_iptm"))
        ptm = c.get("ptm")
        ranking = c.get("confidence_score", c.get("ranking_score"))
        iplddt = norm100(c.get("complex_iplddt", c.get("complex_plddt")))
        # 界面 pAE：优先从 pae npz 精确计算跨链块，否则用 complex_ipde 代理
        pae_npz = find_pae_npz(out_dir, stem)
        if pae_npz:
            i_pae = interface_pae_from_npz(pae_npz, len(target_seq))
        else:
            i_pae = c.get("complex_ipde", "NA")
            if isinstance(i_pae, (int, float)):
                i_pae = round(i_pae, 2)

        def ok(v, thr, ge=True):
            if not isinstance(v, (int, float)):
                return False
            return v >= thr if ge else v <= thr

        passed = int(ok(iptm, args.min_iptm) and
                     (i_pae == "NA" or ok(i_pae, args.max_ipae, ge=False)) and
                     (iplddt == "NA" or ok(iplddt, args.min_iplddt)))
        rows.append({
            "id": sid,
            "iptm": round(iptm, 4) if isinstance(iptm, (int, float)) else iptm,
            "ptm": round(ptm, 4) if isinstance(ptm, (int, float)) else ptm,
            "ranking": round(ranking, 4) if isinstance(ranking, (int, float)) else ranking,
            "interface_pae": i_pae,
            "interface_plddt": iplddt,
            "buried_area_A2": "NA(需SASA工具)",
            "pass": passed,
        })

    # 诊断：若全部为 NA，说明 Boltz 没产出可解析置信度文件，打印目录树辅助排查
    if rows and all(r["iptm"] == "NA" for r in rows):
        print("[警告] 所有候选的 ipTM 均为 NA —— Boltz 可能未产出置信度 JSON。"
              f"输出目录: {out_dir}")
        for root, dirs, files in os.walk(out_dir):
            depth = root[len(out_dir):].count(os.sep)
            if depth > 4:
                continue
            print("  " * depth + os.path.basename(root) + "/")
            for f in sorted(files)[:8]:
                print("  " * (depth + 1) + f)

    # 按 iptm 降序输出
    def sortkey(r):
        return r["iptm"] if isinstance(r["iptm"], (int, float)) else -1
    rows.sort(key=sortkey, reverse=True)

    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "iptm", "ptm", "ranking",
                                           "interface_pae", "interface_plddt",
                                           "buried_area_A2", "pass"])
        w.writeheader()
        w.writerows(rows)

    n_pass = sum(r["pass"] for r in rows)
    print(f"[成功] Boltz 精排 -> {args.out}  (通过 {n_pass}/{len(rows)})")
    print(f"下一步: python 05_rank_candidates.py（融合 ΔG / ipTM 排序）")


if __name__ == "__main__":
    main()
