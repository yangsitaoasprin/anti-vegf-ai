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
02b_rfdiffusion_complex.py
============================================================
抗 VEGF-A 微蛋白 · O3 复合物设计（RFdiffusion partial diffusion + hotspot）

核心思路
--------
直接以 VEGF-A（1FLT 链 W）为靶点，用 RFdiffusion 在其受体结合表面
"原位"扩散生成一条新的 binder 链。输出 PDB 已经是「靶点(A) + binder(B)」
处于结合构象的复合物，可直接喂给 O4 做 MM-GBSA —— 无需单独对接。

这正好解决了先前「单独折叠 binder 硬拼靶点 → ΔG 无物理意义」的问题：
RFdiffusion 在扩散时已知靶点坐标，生成的 binder 天然贴合界面。

用法
----
  # 真实运行（需 RFdiffusion 权重 Complex_base_ckpt.pt）
  python 02b_rfdiffusion_complex.py \
      --target ../data/1FLT.pdb --target-chain W \
      --hotspots ../data/hotspots.txt \
      --lengths 60 80 100 --num-designs 4 \
      --out-dir ../designs/rfdiff

  # 仅验证数据流（不跑模型，生成占位复合物 + 索引）
  python 02b_rfdiffusion_complex.py --demo --out-dir ../designs/rfdiff

依赖: RFdiffusion (D:/RFdiffusion), conda env `esm` (PyTorch CUDA)
"""
import argparse
import csv
import glob
import os
import subprocess
import sys
import tempfile
import site
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.normpath(os.path.join(HERE, ".."))
RFDIFF_DIR = os.environ.get("RFDIFF_DIR", r"D:\RFdiffusion")
ESM_PYTHON = os.environ.get("RFDIFF_PYTHON", r"C:\Users\Administrator\.conda\envs\esm\python.exe")


# ---------------------------------------------------------------------------
# 1. 提取靶点单链（VEGF-A）作为 RFdiffusion 输入
# ---------------------------------------------------------------------------
def extract_chain(pdb_path, chain_id, out_path):
    """提取指定链（保留原始残基编号）写入新 PDB。"""
    n = 0
    with open(pdb_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")) and len(line) > 21 and line[21] == chain_id:
                fout.write(line)
                n += 1
            elif line.startswith(("TER", "END")):
                fout.write(line)
    return n


def chain_residue_set(pdb_path, chain_id):
    s = set()
    for line in open(pdb_path):
        if line.startswith("ATOM") and len(line) > 21 and line[21] == chain_id:
            try:
                s.add(int(line[22:26]))
            except ValueError:
                pass
    return s


def target_range(pdb_path, chain_id):
    res = sorted(chain_residue_set(pdb_path, chain_id))
    if not res:
        return None, None
    return min(res), max(res)


# ---------------------------------------------------------------------------
# 2. 解析 hotspot 文件（过滤水残基与不在靶点链范围内的条目）
# ---------------------------------------------------------------------------
def parse_hotspots(path, target_chain, target_residues):
    """返回 ['W17','W21',...] 形式的 hotspot 列表（仅真实残基）。"""
    out = []
    if not path or not os.path.exists(path):
        return out
    for raw in open(path):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        # 支持两种格式:
        #   "W17  PHE  3.55"   或   "W17"
        tok = line.split()
        if not tok:
            continue
        resid = tok[0]
        # 必须形如 <chain><num>
        if not (resid and resid[0].isalpha() and resid[1:].isdigit()):
            continue
        ch = resid[0]
        num = int(resid[1:])
        # 仅保留靶点链上的真实残基（排除 HOH 等水/配体残基）
        if ch == target_chain and num in target_residues:
            out.append(f"{target_chain}{num}")
    # 去重保序
    seen, uniq = set(), []
    for h in out:
        if h not in seen:
            seen.add(h); uniq.append(h)
    return uniq


# ---------------------------------------------------------------------------
# 3. 精简子进程环境（绕过 Windows 环境变量 >32767 字符导致 CreateProcess 失败）
# ---------------------------------------------------------------------------
def build_env():
    esm = os.path.dirname(ESM_PYTHON)
    lib_bin = os.path.join(esm, "Library", "bin")
    scripts = os.path.join(esm, "Scripts")
    mingw = os.path.join(esm, "Library", "mingw-w64", "bin")
    sys32 = r"C:\Windows\System32"
    sysroot = r"C:\Windows"
    path_parts = [p for p in [scripts, lib_bin, mingw, sys32, sysroot] if os.path.isdir(p)]
    # 追加 CUDA（若存在），torch 自带运行时，但驱动需系统路径
    cuda = r"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.8"
    if os.path.isdir(cuda):
        path_parts.append(os.path.join(cuda, "bin"))
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
        "CUDA_VISIBLE_DEVICES": "0",
        # 解决 Windows 上 libomp.dll 与 libiomp5md.dll 双 OpenMP 运行时冲突
        "KMP_DUPLICATE_LIB_OK": "TRUE",
    }
    # 关键: 必须把用户 site-packages 与 env site-packages 加回 PYTHONPATH,
    # 否则 pandas 等 pip --user 安装的包在精简环境中无法导入
    # (pandas 装在 C:\Users\Administrator\AppData\Roaming\Python\Python311\site-packages)
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


def run_rfdiffusion(rfdiff_dir, input_pdb, out_prefix, length, hotspots,
                    target_chain, rmin, rmax, num_designs, env, cpu=False):
    """调用 RFdiffusion 做复合物（binder）设计。"""
    # 注意: 通过 Python 子进程列表传参时【不要】加外层单引号
    # (官方 shell 示例的单引号仅用于 shell 分词, Hydra 直接解析 [...] 列表)
    contigs = f"[{target_chain}{rmin}-{rmax}/0 {length}-{length}]"
    hot = "[" + ",".join(hotspots) + "]" if hotspots else "null"
    cmd = [
        ESM_PYTHON, os.path.join(rfdiff_dir, "scripts", "run_inference.py"),
        f"inference.input_pdb={input_pdb}",
        f"contigmap.contigs={contigs}",
        f"ppi.hotspot_res={hot}",
        f"inference.num_designs={num_designs}",
        f"inference.output_prefix={out_prefix}",
        "denoiser.noise_scale_ca=0",
        "denoiser.noise_scale_frame=0",
    ]
    print(f"  [RFdiffusion] cmd:\n    {' '.join(cmd)}", flush=True)
    if cpu:
        env = dict(env); env["CUDA_VISIBLE_DEVICES"] = "-1"  # 强制 CPU（DGL Windows 仅 CPU 构建可用）
    r = subprocess.run(cmd, env=env, cwd=rfdiff_dir,
                       capture_output=True, text=True, timeout=3600)
    if r.returncode != 0:
        print(f"  [RFdiffusion][ERROR] rc={r.returncode}")
        print("  ----- stdout -----\n" + r.stdout[-1500:])
        print("  ----- stderr -----\n" + r.stderr[-2000:])
        raise RuntimeError("RFdiffusion 推理失败")
    # 打印关键信息
    for line in r.stdout.splitlines():
        if "Reading models from" in line or "Found GPU" in line or "Finished" in line:
            print("  " + line)
    return r


# ---------------------------------------------------------------------------
# 4. 把 RFdiffusion 输出复合物统一重命名为 A=靶点 / B=binder
# ---------------------------------------------------------------------------
def reassign_complex_chains(pdb_in, pdb_out, target_residues):
    """识别靶点链（残基集合匹配输入靶点）与 binder 链，重命名为 A/B。"""
    per_chain = defaultdict(list)
    for line in open(pdb_in):
        if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
            per_chain[line[21]].append(line)
    # 找出靶点链：其残基集合是输入靶点残基集合的超集
    target_chain = None
    best_overlap = -1
    for c, lines in per_chain.items():
        res = set(int(l[22:26]) for l in lines if len(l) > 26)
        overlap = len(res & target_residues)
        if overlap > best_overlap:
            best_overlap = overlap
            target_chain = c
    if target_chain is None or best_overlap <= 0:
        raise RuntimeError(f"无法在 {pdb_in} 中识别靶点链")
    binder_chain = next((c for c in per_chain if c != target_chain), None)
    if binder_chain is None:
        raise RuntimeError(f"{pdb_in} 中未找到生成的 binder 链")

    written = 0
    with open(pdb_out, "w") as fout:
        for c, lines in per_chain.items():
            new_id = "A" if c == target_chain else "B"
            for l in lines:
                fout.write(l[:21] + new_id + l[22:])
                written += 1
        fout.write("END\n")
    return target_chain, binder_chain, written


# ---------------------------------------------------------------------------
# 5. 主流程
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="RFdiffusion 抗 VEGF-A 复合物(binder)设计")
    ap.add_argument("--target", default=os.path.join(PROJECT, "data", "1FLT.pdb"))
    ap.add_argument("--target-chain", default="W")
    ap.add_argument("--hotspots", default=os.path.join(PROJECT, "data", "hotspots.txt"))
    ap.add_argument("--lengths", nargs="+", type=int, default=[60, 80, 100])
    ap.add_argument("--num-designs", type=int, default=4)
    ap.add_argument("--out-dir", default=os.path.join(PROJECT, "designs", "rfdiff"))
    ap.add_argument("--rfdiffusion-dir", default=RFDIFF_DIR)
    ap.add_argument("--ckpt-dir", default=os.path.join(RFDIFF_DIR, "models"))
    ap.add_argument("--cpu", action="store_true",
                    help="强制 CPU 推理（DGL Windows 仅 CPU 构建可用；慢但稳妥）")
    ap.add_argument("--demo", action="store_true", help="不跑模型，生成占位复合物+索引")
    args = ap.parse_args()

    # 必须转绝对路径: RFdiffusion 以 D:/RFdiffusion 为 cwd 运行,
    # 相对路径会被错误解析（导致 FileNotFoundError）
    args.out_dir = os.path.abspath(args.out_dir)
    os.makedirs(args.out_dir, exist_ok=True)
    target_chain = args.target_chain

    # 解析靶点范围与 hotspots
    rmin, rmax = target_range(args.target, target_chain)
    if rmin is None:
        print(f"[错误] 靶点 {args.target} 中未找到链 {target_chain}")
        sys.exit(1)
    target_residues = chain_residue_set(args.target, target_chain)
    hotspots = parse_hotspots(args.hotspots, target_chain, target_residues)
    print(f"[信息] 靶点链 {target_chain}: 残基 {rmin}-{rmax} (共 {len(target_residues)} 个)")
    print(f"[信息] 有效 hotspot 残基 {len(hotspots)} 个: {hotspots}")

    # 准备单体靶点输入
    target_mono = os.path.join(args.out_dir, f"vegf_A_chain{target_chain}_monomer.pdb")
    n = extract_chain(args.target, target_chain, target_mono)
    print(f"[信息] 提取靶点单链 -> {target_mono} ({n} 个原子行)")

    index_rows = []
    env = build_env()

    for length in args.lengths:
        print(f"\n===== 设计 binder 长度 {length} aa (x{args.num_designs}) =====")
        prefix = os.path.join(args.out_dir, f"vegf_len{length}")
        generated = sorted(glob.glob(f"{prefix}_*.pdb"))

        if args.demo:
            # 占位：构造一个极小复合物（靶点 + 一段随机 binder 坐标）用于验证数据流
            from random import randint
            for i in range(args.num_designs):
                sid = f"vegf_len{length}_{i:03d}"
                out_pdb = os.path.join(args.out_dir, f"{sid}.pdb")
                _write_demo_complex(target_mono, out_pdb, length, randint(100, 999))
                index_rows.append({
                    "id": sid, "length": length, "pdb": os.path.relpath(out_pdb, PROJECT),
                    "hotspot_n": len(hotspots), "note": "[DEMO] 占位复合物",
                })
            continue

        # 真实运行
        if not os.path.exists(os.path.join(args.ckpt_dir, "Complex_base_ckpt.pt")):
            print(f"  [错误] 未找到权重 {os.path.join(args.ckpt_dir, 'Complex_base_ckpt.pt')}")
            print("         请先下载 RFdiffusion 权重 (Complex_base_ckpt.pt)。")
            sys.exit(1)

        if not generated:
            run_rfdiffusion(args.rfdiffusion_dir, target_mono, prefix, length,
                            hotspots, target_chain, rmin, rmax,
                            args.num_designs, env, cpu=args.cpu)
            generated = sorted(glob.glob(f"{prefix}_*.pdb"))

        # 后处理：重命名为 A(靶点)/B(binder)
        for src in generated:
            sid = os.path.splitext(os.path.basename(src))[0]
            dst = os.path.join(args.out_dir, f"{sid}.complex.pdb")
            try:
                tc, bc, w = reassign_complex_chains(src, dst, target_residues)
                index_rows.append({
                    "id": sid, "length": length,
                    "pdb": os.path.relpath(dst, PROJECT),
                    "hotspot_n": len(hotspots),
                    "target_chain": tc, "binder_chain": bc,
                    "note": "RFdiffusion 复合物设计",
                })
                print(f"  [OK] {sid}: 靶点链={tc} binder链={bc} 原子={w}")
            except Exception as e:
                print(f"  [跳过] {sid}: {e}")

    # 写索引
    idx_path = os.path.join(args.out_dir, "index.csv")
    with open(idx_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["id", "length", "pdb", "hotspot_n",
                                          "target_chain", "binder_chain", "note"])
        w.writeheader(); w.writerows(index_rows)
    print(f"\n[完成] 设计索引 -> {idx_path}  (共 {len(index_rows)} 个复合物)")
    print(f"[下一步] 直接喂给 O4（无需对接）:\n"
          f"  python 04_schrodinger_energy.py --complex-dir {os.path.relpath(args.out_dir, HERE)} "
          f"--target-chain A --ligand-chain B --out ../results/energy_real.csv")


def _write_demo_complex(target_mono, out_pdb, length, seed):
    """生成占位复合物（靶点=A 原样 + 一段沿 z 轴排布的 binder=B），仅供数据流测试。"""
    import random
    random.seed(seed)
    with open(out_pdb, "w") as f:
        # 靶点：原样但链改为 A
        for line in open(target_mono):
            if line.startswith(("ATOM", "HETATM")):
                f.write(line[:21] + "A" + line[22:])
            elif line.startswith(("TER", "END")):
                f.write(line)
        # binder：沿 z 轴生成 length 个 CA，链 B
        for i in range(1, length + 1):
            x, y, z = 20.0, 20.0, 5.0 * i
            f.write(
                f"ATOM  {i:5d} CA  ALA B {i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
            )
        f.write("END\n")


if __name__ == "__main__":
    main()
