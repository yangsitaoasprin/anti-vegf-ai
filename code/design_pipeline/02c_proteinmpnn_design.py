# =============================================================================
# Author : Yang Sitao (杨四涛)
# Affiliation : School of Pharmaceutical Sciences, Dali University;
#               The Third Affiliated Hospital of Dali University
# Contact : yangsitaoasprin@126.com
# License : MIT (see ../../LICENSE)
# =============================================================================
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02c_proteinmpnn_design.py
============================================================
抗 VEGF-A 微蛋白 · 序列设计（ProteinMPNN）

作用
----
把 RFdiffusion 复合物设计（02b）产出的骨架复合物（A=靶点 VEGF-A、B=binder）
喂给 ProteinMPNN，**固定靶点链 A、只为 binder 链 B 设计氨基酸序列**，
多温度采样产生多样化候选序列，输出 FASTA 供下游 ESMFold 粗筛（03）。

流程（全部用精简环境 subprocess 调 ProteinMPNN，绕过 Windows 环境变量超长坑）
----
  1) 收集 complex-dir 下的 *.complex.pdb 到工作目录
  2) parse_multiple_chains.py  → parsed_chains.jsonl
  3) assign_fixed_chains.py --chain_list "B"  → 固定 A、设计 B
  4) protein_mpnn_run.py 多温度采样
  5) 提取每条设计的 **binder 链(B)** 序列，写成单序列 FASTA 到 out-dir

用法
----
  python 02c_proteinmpnn_design.py \
      --complex-dir ../designs/rfdiff_gpu \
      --out-dir ../designs/mpnn \
      --num-seq 8 --temps 0.1 0.2 0.3

  # 仅验证数据流（不跑模型）
  python 02c_proteinmpnn_design.py --complex-dir ../designs/rfdiff_gpu --demo

依赖: D:/ProteinMPNN（含 vanilla_model_weights），conda env `esm`（PyTorch，GPU/CPU 均可）
"""
import argparse
import csv
import glob
import os
import shutil
import site
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.normpath(os.path.join(HERE, ".."))
MPNN_DIR = os.environ.get("MPNN_DIR", r"D:\ProteinMPNN")
# ProteinMPNN 只依赖 torch（不需 DGL），esm 环境的 GPU torch 即可
MPNN_PYTHON = os.environ.get("MPNN_PYTHON",
                             r"C:\Users\Administrator\.conda\envs\esm\python.exe")


# ---------------------------------------------------------------------------
# 精简子进程环境（绕过 Windows 环境变量 >32767 字符导致 CreateProcess 失败）
# ---------------------------------------------------------------------------
def build_env():
    envdir = os.path.dirname(MPNN_PYTHON)
    lib_bin = os.path.join(envdir, "Library", "bin")
    scripts = os.path.join(envdir, "Scripts")
    mingw = os.path.join(envdir, "Library", "mingw-w64", "bin")
    sys32 = r"C:\Windows\System32"
    sysroot = r"C:\Windows"
    path_parts = [p for p in [envdir, scripts, lib_bin, mingw, sys32, sysroot]
                  if os.path.isdir(p)]
    env = {
        "PATH": ";".join(path_parts),
        "SYSTEMROOT": sysroot,
        "WINDIR": sysroot,
        "SYSTEMDRIVE": "C:",
        "COMSPEC": os.path.join(sys32, "cmd.exe"),
        "TEMP": os.environ.get("TEMP", r"C:\Windows\Temp"),
        "TMP": os.environ.get("TMP", r"C:\Windows\Temp"),
        "USERPROFILE": os.environ.get("USERPROFILE", r"C:\Users\Administrator"),
        "HOMEDRIVE": "C:",
        "HOMEPATH": r"\Users\Administrator",
        "USERNAME": os.environ.get("USERNAME", "Administrator"),
        "OS": "Windows_NT",
        "PROCESSOR_ARCHITECTURE": "AMD64",
        "NUMBER_OF_PROCESSORS": os.environ.get("NUMBER_OF_PROCESSORS", "8"),
        "PYTHONPATH": "",
        "PYTHONIOENCODING": "utf-8",
        "CUDA_VISIBLE_DEVICES": os.environ.get("CUDA_VISIBLE_DEVICES", "0"),
        "KMP_DUPLICATE_LIB_OK": "TRUE",
    }
    # 把 site-packages 加回 PYTHONPATH（pip --user 装的包在精简环境里也能导入）
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
    return env


def run(cmd, env, cwd=None):
    print("  $ " + " ".join(str(c) for c in cmd), flush=True)
    r = subprocess.run(cmd, env=env, cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        print("  [stderr]\n" + (r.stderr or "")[-1500:], flush=True)
    return r


# ---------------------------------------------------------------------------
# FASTA 解析：从 ProteinMPNN 输出提取 binder(B) 链序列
# ---------------------------------------------------------------------------
def parse_mpnn_fasta(fa_path, binder_index=1):
    """
    ProteinMPNN 输出 FASTA（固定链A、设计链B 时）：
      第一条记录 = 原生序列（binder 链B的脚手架占位，跳过）
      后续每条记录 = 单个设计链B的完整序列（本机构型只输出被设计的链）
    因此直接取每条非原生记录的完整序列作为 binder 序列即可，无需按 '/' 切分。
    返回 [(tag, binder_seq), ...]（跳过第一条原生序列）。
    """
    recs = []
    header, seq = None, []
    with open(fa_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    recs.append((header, "".join(seq)))
                header, seq = line[1:], []
            elif line:
                seq.append(line.strip())
    if header is not None:
        recs.append((header, "".join(seq)))
    AA = set("ACDEFGHIKLMNPQRSTVWY")
    entries = []
    for idx, (hdr, full) in enumerate(recs):
        if idx == 0:
            continue  # 原生序列
        s = "".join(c for c in full.upper() if c in AA)
        if len(s) < 10:
            continue
        entries.append((hdr, s))
    return entries


def main():
    ap = argparse.ArgumentParser(description="ProteinMPNN 序列设计（固定靶点A、设计binderB）")
    ap.add_argument("--complex-dir", required=True,
                    help="RFdiffusion 复合物目录（含 *.complex.pdb，A=靶点/B=binder）")
    ap.add_argument("--out-dir", default=os.path.join(PROJECT, "designs", "mpnn"),
                    help="输出 FASTA 目录")
    ap.add_argument("--design-chain", default="B", help="设计链（binder），默认 B")
    ap.add_argument("--num-seq", type=int, default=8, help="每骨架每温度序列数")
    ap.add_argument("--temps", nargs="+", default=["0.1", "0.2", "0.3"],
                    help="采样温度列表")
    ap.add_argument("--model-name", default="v_48_020",
                    help="ProteinMPNN 模型：v_48_002/010/020/030")
    ap.add_argument("--seed", type=int, default=37)
    ap.add_argument("--demo", action="store_true", help="不跑模型，仅验证数据流")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    complexes = sorted(glob.glob(os.path.join(args.complex_dir, "*.complex.pdb")))
    print(f"[信息] 复合物骨架 {len(complexes)} 个 (来自 {args.complex_dir})")
    if not complexes:
        print("[错误] 未找到 *.complex.pdb，请先运行 02b_rfdiffusion_complex.py")
        sys.exit(1)

    if args.demo:
        print("[DEMO] 数据流验证（不跑 ProteinMPNN）:")
        for c in complexes[:5]:
            sid = os.path.basename(c).replace(".complex.pdb", "")
            print(f"  {sid}: 固定链 A(靶点) / 设计链 {args.design_chain}(binder), "
                  f"{args.num_seq}×{len(args.temps)} 条序列")
        print(f"[DEMO] 将输出到 {args.out_dir}/<id>.fa")
        return

    if not os.path.isdir(MPNN_DIR):
        print(f"[错误] 未找到 ProteinMPNN: {MPNN_DIR}（设 MPNN_DIR 环境变量指向安装目录）")
        sys.exit(1)

    env = build_env()
    weights = os.path.join(MPNN_DIR, "vanilla_model_weights")
    parse_py = os.path.join(MPNN_DIR, "helper_scripts", "parse_multiple_chains.py")
    assign_py = os.path.join(MPNN_DIR, "helper_scripts", "assign_fixed_chains.py")
    run_py = os.path.join(MPNN_DIR, "protein_mpnn_run.py")

    work = tempfile.mkdtemp(prefix="mpnn_")
    pdb_in = os.path.join(work, "pdbs")
    os.makedirs(pdb_in, exist_ok=True)
    for c in complexes:
        shutil.copy(c, os.path.join(pdb_in, os.path.basename(c).replace(".complex", "")))

    parsed = os.path.join(work, "parsed_chains.jsonl")
    chain_id = os.path.join(work, "assigned_chains.jsonl")
    mpnn_out = os.path.join(work, "mpnn_out")
    os.makedirs(mpnn_out, exist_ok=True)

    print("[1/3] 解析复合物链结构 ...")
    r = run([MPNN_PYTHON, parse_py, f"--input_path={pdb_in}",
             f"--output_path={parsed}"], env)
    if r.returncode != 0 or not os.path.exists(parsed):
        print("[错误] parse_multiple_chains 失败"); sys.exit(2)

    print(f"[2/3] 指定设计链={args.design_chain}（固定靶点链 A）...")
    r = run([MPNN_PYTHON, assign_py, f"--input_path={parsed}",
             f"--output_path={chain_id}", f"--chain_list={args.design_chain}"], env)
    if r.returncode != 0 or not os.path.exists(chain_id):
        print("[错误] assign_fixed_chains 失败"); sys.exit(3)

    print(f"[3/3] ProteinMPNN 采样：温度={args.temps} 每骨架每温度 {args.num_seq} 条 ...")
    # 注意：本机 ProteinMPNN 的 protein_mpnn_run.py 的 --sampling_temp 仅接受
    # 单个温度，故按温度分别跑一次，各自输出到独立子目录避免文件名冲突。
    fa_files = []
    for t in args.temps:
        t_out = os.path.join(work, f"mpnn_out_{t}")
        os.makedirs(t_out, exist_ok=True)
        r = run([MPNN_PYTHON, run_py,
                 "--jsonl_path", parsed,
                 "--chain_id_jsonl", chain_id,
                 "--out_folder", t_out,
                 "--num_seq_per_target", str(args.num_seq),
                 "--sampling_temp", t,
                 "--batch_size", "1",
                 "--seed", str(args.seed),
                 "--path_to_model_weights", weights,
                 "--model_name", args.model_name], env)
        if r.returncode != 0:
            print(f"[错误] protein_mpnn_run 失败 (temp={t})"); sys.exit(4)
        seq_dir = os.path.join(t_out, "seqs")
        if os.path.isdir(seq_dir):
            # 记录 (路径, 温度)，避免不同温度产物同名覆盖
            fa_files.extend((f, t) for f in sorted(glob.glob(os.path.join(seq_dir, "*.fa"))))
    print(f"[信息] ProteinMPNN 产出 {len(fa_files)} 个 FASTA，提取 binder 链序列 ...")

    index_rows = []
    total = 0
    for fa, t in fa_files:
        base = os.path.splitext(os.path.basename(fa))[0]  # e.g. vegf_len60_0
        entries = parse_mpnn_fasta(fa)
        for k, (hdr, binder_seq) in enumerate(entries):
            sid = f"{base}_T{t}_s{k}"
            out_fa = os.path.join(args.out_dir, f"{sid}.fa")
            with open(out_fa, "w", encoding="utf-8") as w:
                w.write(f">{sid} src={base} {hdr}\n{binder_seq}\n")
            index_rows.append({"id": sid, "source_backbone": base,
                               "length": len(binder_seq), "mpnn_header": hdr})
            total += 1

    idx_csv = os.path.join(args.out_dir, "mpnn_index.csv")
    with open(idx_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "source_backbone", "length", "mpnn_header"])
        w.writeheader()
        w.writerows(index_rows)

    try:
        shutil.rmtree(work)
    except Exception:
        pass

    print(f"[成功] 共 {total} 条 binder 序列 -> {args.out_dir}")
    print(f"       索引: {idx_csv}")
    print(f"下一步: python 03_esmfold_prescreen.py --seqs {args.out_dir} "
          f"--model-path D:/pmodel/models/esmfold_v1")


if __name__ == "__main__":
    main()
