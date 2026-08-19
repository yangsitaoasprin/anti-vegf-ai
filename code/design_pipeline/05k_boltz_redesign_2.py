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
05k_boltz_redesign_2.py — #2 去A重设计候选的 Boltz-1 精排
对 05i ESMFold 预筛通过者（results/deA_redesign_1_2/#2_esmfold_passed.fa）做
Boltz-1(AF3替代) 复合物预测，解析 ipTM / ptm / 界面pAE / 界面pLDDT，按与基线
#2 (vegf_len100_4, baseline ipTM≈0.63 / 界面pAE≈8.5) 同口径排名。

靶点 = #2 设计靶点（complex PDB chain A，98aa VEGF-A 受体结合核心），
与 boltz_rerank.csv 基线完全同口径 → 结果可直接比较。

关键补丁（继承 04）:
  --model boltz1 --no_kernels  (RTX5090 sm_120 必加，否则 cuequivariance 算子缺失)

用法: python 05k_boltz_redesign_2.py [--seeds 0] [--limit N] [--out CSV]
"""
import os, sys, subprocess, glob, csv, shutil, tempfile, re

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, ".."))
BOLTZ_PYTHON = r"C:\Users\Administrator\.conda\envs\boltz_cuda\python.exe"
BOLTZ_CLI = os.path.join(os.path.dirname(BOLTZ_PYTHON), "Scripts", "boltz.exe")

# #2 设计靶点 = 原始 complex PDB chain A（98aa VEGF-A RBD 核心），与基线完全一致
TARGET = ("HEVVKFMDVYQRSYCHPIETLVDIFQEYPDEIEYIFKPSCVPLMRCGGCCNDEGLECVPTEESNIT"
          "MQIMRIKPHQGQHIGEMSFLQHNKCECRPKKD")

IN_FA = os.path.join(PROJ, "results", "deA_redesign_1_2", "#2_esmfold_passed.fa")
OUT_CSV = os.path.join(PROJ, "results", "deA_redesign_1_2", "boltz_redesign_2.csv")
WORK_ROOT = os.path.join(PROJ, "results", "deA_redesign_1_2", "boltz_pdb")
LEN_A = len(TARGET)

# 基线对比（来自 results/boltz_rerank.csv 的 len100_4 最优）
BASE_IPTM = 0.63
BASE_IPAE = 8.5
BASE_IPLDDT = 77.7


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


def write_yaml(path, binder_seq, seed):
    with open(path, "w", encoding="utf-8") as f:
        f.write("version: 1\n")
        f.write("sequences:\n")
        f.write("  - protein:\n")
        f.write("      id: A\n")
        f.write(f"      sequence: {TARGET}\n")
        f.write("      msa: empty\n")
        f.write("  - protein:\n")
        f.write("      id: B\n")
        f.write(f"      sequence: {binder_seq}\n")
        f.write("      msa: empty\n")


def find_confidence_json(out_dir, stem):
    # 注意：stem 内含 [K95Q,A85L] 等方括号，glob 会把 [] 当字符类误判，
    # 故用 os.walk 精确前缀匹配，避开 glob 的字符类陷阱。
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
    hits = glob.glob(os.path.join(out_dir, "**", base), recursive=True)
    return sorted(hits)[0] if hits else None


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


def read_fastas(path):
    recs, name = [], None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
                recs.append([name, ""])
            elif name:
                recs[-1][1] += line
    return [(n, s) for n, s in recs if s]


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_fa", default=IN_FA)
    ap.add_argument("--out", default=OUT_CSV)
    ap.add_argument("--seeds", default="0", help="逗号分隔的随机种子，如 0,11,42")
    ap.add_argument("--limit", type=int, default=0, help="最多处理前 N 条")
    args = ap.parse_args()

    seeds = [int(x) for x in args.seeds.split(",") if x != ""]
    seqs = read_fastas(args.in_fa)
    if args.limit:
        seqs = seqs[:args.limit]
    print(f"[信息] 待 Boltz 精排: {len(seqs)} 条 (seeds={seeds}) | 靶点A={LEN_A}aa")

    env = build_env()
    os.makedirs(WORK_ROOT, exist_ok=True)
    rows = []

    for si, (name, seq) in enumerate(seqs, 1):
        seed_metrics = []  # (seed, iptm, ptm, ipae, iplddt)
        for seed in seeds:
            stem = f"{name}_s{seed}"
            in_dir = os.path.join(WORK_ROOT, stem, "inputs")
            out_dir = os.path.join(WORK_ROOT, stem, "out")
            os.makedirs(in_dir, exist_ok=True)
            # 同 candidate 多次 seed 用不同 stem，避免覆盖
            yaml_stem = stem
            yaml_path = os.path.join(in_dir, f"{yaml_stem}.yaml")
            write_yaml(yaml_path, seq, seed)

            cmd = [BOLTZ_CLI, "predict", in_dir, "--out_dir", out_dir,
                   "--output_format", "pdb", "--recycling_steps", "3",
                   "--write_full_pae", "--model", "boltz1", "--no_kernels",
                   "--override", "--seed", str(seed)]
            print(f"[{si}/{len(seqs)}] {name} seed={seed}  $ boltz predict ...", flush=True)
            r = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=1200)
            if r.returncode != 0:
                print(f"  FAIL seed={seed}: {r.stderr[-600:]}", flush=True)
                seed_metrics.append((seed, None, None, None, None))
                continue
            cj = find_confidence_json(out_dir, yaml_stem)
            iptm = ptm = iplddt = None
            if cj:
                import json
                d = json.load(open(cj))
                # iptm/ptm 已是 0-1 尺度，切勿 ×100（否则 0.13→13.01 伪值）
                iptm = d.get("iptm")
                ptm = d.get("ptm")
                # 仅 pLDDT 类需归一到 0-100
                iplddt = d.get("complex_iplddt") if d.get("complex_iplddt") is not None \
                         else d.get("complex_plddt")
            ipae = "NA"
            npz = find_pae_npz(out_dir, yaml_stem)
            if npz:
                ipae = interface_pae_from_npz(npz, LEN_A)
            # 复制 PDB 到扁平目录
            pdbs = glob.glob(os.path.join(out_dir, "**", f"{yaml_stem}_model_0.pdb"), recursive=True)
            if pdbs:
                shutil.copy2(pdbs[0], os.path.join(WORK_ROOT, os.path.basename(pdbs[0])))
            # 仅 pLDDT 需归一到 0-100；iptm/ptm 保持 0-1
            iplddt = norm100(iplddt) if iplddt is not None else None
            seed_metrics.append((seed, iptm, ptm, ipae, iplddt))
            print(f"  -> seed={seed} ipTM={iptm} ptm={ptm} 界面pAE={ipae} 界面pLDDT={iplddt}", flush=True)

        # 多 seed 取均值（ipTM/界面pAE/界面pLDDT）
        def avg(vals):
            vs = [v for v in vals if isinstance(v, (int, float))]
            return round(sum(vs) / len(vs), 3) if vs else None
        iptms = [m[1] for m in seed_metrics]
        paes = [m[3] for m in seed_metrics if isinstance(m[3], (int, float))]
        pls = [m[4] for m in seed_metrics if isinstance(m[4], (int, float))]
        m_iptm = avg(iptms)
        m_pae = avg(paes)
        m_pl = avg(pls)
        rows.append({
            "id": name, "len": len(seq),
            "iptm_mean": m_iptm if m_iptm is not None else "NA",
            "iptm_best": max([v for v in iptms if isinstance(v, (int, float))], default="NA"),
            "ptm_mean": avg([m[2] for m in seed_metrics]),
            "interface_pae_mean": m_pae if m_pae is not None else "NA",
            "interface_plddt_mean": m_pl if m_pl is not None else "NA",
            "n_seed": len(seeds),
            "pass": int(m_iptm >= 0.5) if m_iptm is not None else 0,
            "vs_baseline": f"基线ipTM{BASE_IPTM}/pAE{BASE_IPAE}",
        })

    # 排序：ipTM 降序，界面pAE 升序
    def sortkey(r):
        ip = r["iptm_mean"]
        pa = r["interface_pae_mean"]
        ip = ip if isinstance(ip, (int, float)) else -1
        pa = pa if isinstance(pa, (int, float)) else 999
        return (-ip, pa)
    rows.sort(key=sortkey)

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["id", "len", "iptm_mean", "iptm_best",
                              "ptm_mean", "interface_pae_mean", "interface_plddt_mean",
                              "n_seed", "pass", "vs_baseline"])
        w.writeheader()
        w.writerows(rows)
    print(f"\n[成功] Boltz 精排 -> {args.out} ({len(rows)} 条)")
    print("=== 排名（ipTM降序 / 界面pAE升序）===")
    for r in rows:
        print(f"  {r['id']:20s} ipTM={r['iptm_mean']} ptm={r['ptm_mean']} "
              f"界面pAE={r['interface_pae_mean']} 界面pLDDT={r['interface_plddt_mean']} "
              f"pass={r['pass']}")
    print(f"\n基线 #2 最优: ipTM={BASE_IPTM} 界面pAE={BASE_IPAE} 界面pLDDT={BASE_IPLDDT}")


if __name__ == "__main__":
    main()
