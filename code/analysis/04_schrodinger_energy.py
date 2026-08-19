# =============================================================================
# Author : Yang Sitao (杨四涛)
# Affiliation : School of Pharmaceutical Sciences, Dali University;
#               The Third Affiliated Hospital of Dali University
# Contact : yangsitaoasprin@126.com
# License : MIT (see ../../LICENSE)
# =============================================================================
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")
"""
04_schrodinger_energy.py
界面能量评估：使用 Schrödinger Prime MM-GBSA 计算候选微蛋白与 VEGF-A 的结合自由能。

正确的蛋白-蛋白工作流：
    1. 准备复合物结构 (prepwizard: 加氢 / 质子化 / 优化 H-bond 网络)
    2. [可选但强烈推荐] 蛋白-蛋白对接 (prime_zbind) 产生「真实结合构象」
    3. prime_mmgbsa 计算 MM-GBSA ΔG_bind
    4. 解析 r_psp_MMGBSA_dG_Bind 及其分量，写入 CSV

⚠ 关键科学说明
    MM-GBSA 必须作用在「真实结合构象」上。若 binder 是单独折叠后
    按原始坐标硬拼到靶点（未经对接），所得 ΔG 无物理意义，典型表现为
    巨大正值（如 +3000 kcal/mol）。
    生产流水线中 O3 (RFdiffusion) 直接在靶点结合位点内设计 binder，
    输出即复合物构象，此时本模块 MM-GBSA 才有效；若输入为单独折叠的
    binder，则应先经 prime_zbind 对接（本脚本 --dock 选项）。

用法：
    # ── 直接喂 RFdiffusion 复合物设计输出（推荐，无需对接）──
    # 复合物已是「靶点(A) + binder(B)」结合构象，prepwizard 后直接 MM-GBSA
    python 04_schrodinger_energy.py --complex-dir ../designs/rfdiff \
        --target-chain A --ligand-chain B --out ../results/energy_real.csv

    # 单独折叠 binder → 先对接再算能量
    python 04_schrodinger_energy.py --pdb-dir ../data/screen_esmfold \
        --target ../data/1FLT.pdb --target-chain W --dock --out ../results/energy.csv

    # demo 模式（无 Schrödinger 时演示数据流）
    python 04_schrodinger_energy.py --demo --out ../results/energy.csv

依赖：Schrödinger 2025-2 (prepwizard, prime_mmgbsa；对接可选 prime_zbind)
"""

import argparse
import csv
import glob
import os
import subprocess
import sys
import tempfile

SCHRODINGER = r"C:\Program Files\Schrodinger2025-2"
STRUCTCONVERT = os.path.join(SCHRODINGER, "utilities", "structconvert.exe")
PREPWIZARD = os.path.join(SCHRODINGER, "utilities", "prepwizard.exe")
PRIME_MMGBSA = os.path.join(SCHRODINGER, "prime_mmgbsa.exe")
PRIME_ZBIND = os.path.join(SCHRODINGER, "prime_zbind.exe")

# ---------------------------------------------------------------------------
# 干净环境启动：绕过 WorkBuddy 会话中 >32767 字符的环境变量导致
# Schrödinger toplevel.py 在 os.environ.update() 崩溃的问题
# ---------------------------------------------------------------------------

def _clean_env():
    return {
        "PATH": (
            r"C:\Windows\System32;C:\Windows;" + SCHRODINGER + ";"
            + os.path.join(SCHRODINGER, "utilities") + ";"
            + os.path.join(SCHRODINGER, "mmshare-v7.0", "bin", "Windows-x64")
        ),
        "SCHRODINGER": SCHRODINGER,
        "SYSTEMROOT": r"C:\Windows",
        "SystemRoot": r"C:\Windows",
        "COMSPEC": r"C:\Windows\System32\cmd.exe",
        "TEMP": os.environ.get("TEMP", r"C:\Windows\Temp"),
        "TMP": os.environ.get("TMP", r"C:\Windows\Temp"),
        "NUMBER_OF_PROCESSORS": os.environ.get("NUMBER_OF_PROCESSORS", "8"),
        "PROCESSOR_ARCHITECTURE": "AMD64",
        "OS": "Windows_NT",
        "USERPROFILE": os.environ.get("USERPROFILE", r"C:\Users\Administrator"),
        "HOMEDRIVE": "C:",
        "HOMEPATH": r"\Users\Administrator",
        "USERNAME": os.environ.get("USERNAME", "Administrator"),
        "WINDIR": r"C:\Windows",
    }


def _run_tool(tool_path, args, cwd):
    """用干净环境运行 Schrödinger 工具，返回 CompletedProcess。"""
    return subprocess.run(
        [tool_path] + args, env=_clean_env(),
        capture_output=True, text=True, cwd=cwd, timeout=1800,
    )


# ---------------------------------------------------------------------------
# 结构处理工具
# ---------------------------------------------------------------------------

def check_schrodinger():
    ok = True
    for name, path in [("structconvert", STRUCTCONVERT), ("prime_mmgbsa", PRIME_MMGBSA)]:
        if not os.path.exists(path):
            print(f"[错误] 未找到 {name}: {path}")
            ok = False
    return ok


def extract_chain(pdb_path, chain_id, out_path):
    """从 PDB 提取指定链到新文件（保留原链列）。"""
    with open(pdb_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")) and len(line) > 21 and line[21] == chain_id:
                fout.write(line)
            elif line.startswith(("TER", "END")):
                fout.write(line)


def merge_pdb(receptor_pdb, ligand_pdb, out_path, receptor_chain="A", ligand_chain="B"):
    """合并受体与配体 PDB 为复合物，重命名链 ID。"""
    with open(out_path, "w") as fout:
        with open(receptor_pdb) as fin:
            for line in fin:
                if line.startswith(("ATOM", "HETATM")):
                    line = line[:21] + receptor_chain + line[22:]
                    fout.write(line)
                elif line.startswith("TER"):
                    fout.write(line)
        with open(ligand_pdb) as fin:
            for line in fin:
                if line.startswith(("ATOM", "HETATM")):
                    line = line[:21] + ligand_chain + line[22:]
                    fout.write(line)
                elif line.startswith("TER"):
                    fout.write(line)
        fout.write("END\n")


def convert_to_mae(pdb_path, mae_path, work_dir):
    r = _run_tool(STRUCTCONVERT, [pdb_path, mae_path], work_dir)
    if r.returncode != 0 or not os.path.exists(mae_path):
        raise RuntimeError(f"structconvert 失败 (rc={r.returncode}): {r.stderr[-500:]}")


def run_prepwizard(mae_in, mae_out, work_dir):
    """prepwizard 准备结构：加氢、质子化、优化。输入可为 PDB 或 MAE。"""
    mae_in = os.path.abspath(mae_in)
    mae_out = os.path.abspath(mae_out)
    r = _run_tool(PREPWIZARD, [
        mae_in, mae_out,
        "-fillsidechains", "-disulfides", "-propka_pH", "7.4",
        "-HOST", "localhost", "-WAIT",
    ], work_dir)
    if r.returncode != 0 or not os.path.exists(mae_out):
        raise RuntimeError(f"prepwizard 失败 (rc={r.returncode}): {r.stderr[-800:]}")
    return r


def run_zbind(mae_in, mae_out, work_dir):
    """prime_zbind 蛋白-蛋白对接，产生结合构象。返回 CompletedProcess。"""
    return _run_tool(PRIME_ZBIND, [
        mae_in, mae_out,
        "-poses", "10", "-HOST", "localhost", "-WAIT",
    ], work_dir)


def run_mmgbsa(mae_path, work_dir, ligand_asl):
    """prime_mmgbsa 计算 MM-GBSA 结合自由能。
    输出 CSV 以输入文件名派生：<stem>-out.csv（位于 work_dir）。"""
    mae_path = os.path.abspath(mae_path)
    r = _run_tool(PRIME_MMGBSA, [
        mae_path,
        "-job_type", "ENERGY",
        "-csv_output", "yes",
        "-ligand", ligand_asl,
        "-HOST", "localhost", "-WAIT",
    ], work_dir)
    return r


def find_mmgbsa_csv(work_dir, stem):
    """定位 prime_mmgbsa 输出 CSV。"""
    cand = os.path.join(work_dir, f"{stem}-out.csv")
    if os.path.exists(cand):
        return cand
    matches = glob.glob(os.path.join(work_dir, "*-out.csv"))
    return matches[0] if matches else None


def parse_mmgbsa_csv(csv_path):
    """解析 MM-GBSA CSV，返回 dG 及分量字典。"""
    out = {}
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        row = None
        for row in reader:
            break
        if row is None:
            return out
        for k, v in row.items():
            if k.startswith("r_psp_MMGBSA_dG_Bind") or k.startswith("r_psp_MMGBSA_dG_Bind("):
                try:
                    out[k] = float(v)
                except (ValueError, TypeError):
                    out[k] = v
    return out


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------

def process_candidate(binder_pdb, target_pdb, target_chain, ligand_chain,
                      work_dir, job_name, do_dock):
    """处理单个候选：合并 → 准备 → [对接] → MM-GBSA。"""
    os.makedirs(work_dir, exist_ok=True)

    # 1. 提取靶点链
    receptor_pdb = os.path.join(work_dir, "receptor.pdb")
    extract_chain(target_pdb, target_chain, receptor_pdb)

    # 2. 合并复合物
    complex_pdb = os.path.join(work_dir, f"{job_name}_complex.pdb")
    merge_pdb(receptor_pdb, binder_pdb, complex_pdb,
              receptor_chain="A", ligand_chain=ligand_chain)

    # 3. 准备结构
    prepared_mae = os.path.join(work_dir, f"{job_name}_prepared.mae")
    run_prepwizard(complex_pdb, prepared_mae, work_dir)

    # 4. 可选：对接产生真实结合构象
    docked = False
    mae_for_mmgbsa = prepared_mae
    if do_dock and os.path.exists(PRIME_ZBIND):
        docked_mae = os.path.join(work_dir, f"{job_name}_docked.mae")
        rz = run_zbind(prepared_mae, docked_mae, work_dir)
        if rz.returncode == 0 and os.path.exists(docked_mae):
            mae_for_mmgbsa = docked_mae
            docked = True
        else:
            print(f"    [提示] prime_zbind 不可用/失败，退回未经对接的复合物（ΔG 可能无物理意义）")

    # 5. MM-GBSA
    ligand_asl = f'chain.name {ligand_chain}'
    rmm = run_mmgbsa(mae_for_mmgbsa, work_dir, ligand_asl)
    stem = os.path.splitext(os.path.basename(mae_for_mmgbsa))[0]
    csv_path = find_mmgbsa_csv(work_dir, stem)
    energies = parse_mmgbsa_csv(csv_path) if csv_path else {}

    return {
        "docked": docked,
        "prepared_mae": prepared_mae,
        "mae_for_mmgbsa": mae_for_mmgbsa,
        "mmgbsa_stdout": rmm.stdout[-800:] if rmm else "",
        "mmgbsa_stderr": rmm.stderr[-800:] if rmm else "",
        "energies": energies,
    }


def count_chain_residues(pdb_path, chain_id):
    """统计某链残基数。"""
    res = set()
    for line in open(pdb_path):
        if line.startswith("ATOM") and len(line) > 21 and line[21] == chain_id:
            try:
                res.add(int(line[22:26]))
            except ValueError:
                pass
    return len(res)


def has_chain(pdb_path, chain_id):
    for line in open(pdb_path):
        if line.startswith(("ATOM", "HETATM")) and len(line) > 21 and line[21] == chain_id:
            return True
    return False


def process_complex(complex_pdb, target_chain, ligand_chain, work_dir, job_name):
    """直接处理 RFdiffusion 复合物（靶点+已就位 binder），免对接。

    输入复合物已是结合构象：靶点链=target_chain，binder链=ligand_chain。
    仅做 prepwizard 准备 → prime_mmgbsa（ligand=ligand_chain）。
    """
    os.makedirs(work_dir, exist_ok=True)

    # 1. 准备复合物结构
    prepared_mae = os.path.join(work_dir, f"{job_name}_prepared.mae")
    run_prepwizard(complex_pdb, prepared_mae, work_dir)

    # 2. MM-GBSA（ligand = binder 链，其余=受体）
    ligand_asl = f'chain.name {ligand_chain}'
    rmm = run_mmgbsa(prepared_mae, work_dir, ligand_asl)
    stem = os.path.splitext(os.path.basename(prepared_mae))[0]
    csv_path = find_mmgbsa_csv(work_dir, stem)
    energies = parse_mmgbsa_csv(csv_path) if csv_path else {}

    return {
        "prepared_mae": prepared_mae,
        "mmgbsa_stdout": rmm.stdout[-800:] if rmm else "",
        "mmgbsa_stderr": rmm.stderr[-800:] if rmm else "",
        "energies": energies,
    }


def main():
    ap = argparse.ArgumentParser(description="Schrödinger Prime MM-GBSA 界面能量评估")
    ap.add_argument("--candidates", default=None, help="粗筛结果 CSV")
    ap.add_argument("--pdb-dir", default=None, help="候选 PDB 目录")
    ap.add_argument("--target", default="../data/1FLT.pdb", help="靶点 PDB")
    ap.add_argument("--target-chain", default="W", help="靶点链 ID")
    ap.add_argument("--ligand-chain", default="B", help="合并后 binder 链 ID（ASL 用）")
    ap.add_argument("--pass-only", action="store_true", help="仅处理通过粗筛的候选")
    ap.add_argument("--top-n", type=int, default=None, help="仅处理前 N 个候选")
    ap.add_argument("--dock", action="store_true", help="先用 prime_zbind 对接产生结合构象")
    ap.add_argument("--out", default="../results/energy.csv")
    ap.add_argument("--demo", action="store_true", help="演示模式")
    ap.add_argument("--complex-dir", default=None,
                    help="直接喂 RFdiffusion 复合物 PDB 目录（免对接）。"
                         "默认处理其中含 --ligand-chain 的 *.pdb")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)

    # Demo 模式
    if args.demo:
        print("[DEMO] 模拟 Schrödinger Prime MM-GBSA 能量评估")
        rows = []
        for i in range(5):
            rows.append({
                "id": f"vegf_binder_{i:03d}",
                "length": 60 + i * 5,
                "mean_plddt": 85.0 - i * 2,
                "prime_dg": round(-8.5 + i * 0.8, 2),
                "dg_hbond": round(-2.1 + i * 0.3, 2),
                "dg_vdw": round(-4.2 + i * 0.4, 2),
                "dg_solv": round(-1.5 + i * 0.2, 2),
                "docked": 0,
                "pass": 1 if i < 3 else 0,
                "note": "[DEMO] 占位值，非真实计算",
            })
        with open(args.out, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"[DEMO] 能量结果 -> {args.out}")
        print(f"[DEMO] 第一名: {rows[0]['id']}  ΔG={rows[0]['prime_dg']} kcal/mol")
        return

    if not check_schrodinger():
        print("\n[指引] 未检测到 Schrödinger。可选项:")
        print("  1. 安装 Schrödinger 2025-2 (商业许可，本项目已有)")
        print("  2. 使用 FoldX (学术免费, 原生 Windows)")
        print("  3. 加 --demo 查看数据流示例")
        sys.exit(1)

    # ── 直接喂 RFdiffusion 复合物（无需对接）──
    if args.complex_dir:
        print(f"[信息] 复合物直算模式（免对接）: {args.complex_dir}")
        print(f"[信息] 靶点链={args.target_chain}  binder链={args.ligand_chain}")
        cands = []
        # 只取 O3 后处理生成的 *.complex.pdb（已重命名为 A=靶点/B=binder），
        # 跳过原始 RFdiffusion 输出(靶点链W+自动命名binder)与靶点单体
        for pdb_path in sorted(glob.glob(os.path.join(args.complex_dir, "*.complex.pdb"))):
            pdb_path = os.path.abspath(pdb_path)
            sid = os.path.splitext(os.path.basename(pdb_path))[0]
            if not has_chain(pdb_path, args.ligand_chain):
                print(f"  [跳过] {sid}: 不含 binder 链 {args.ligand_chain}")
                continue
            cands.append((sid, pdb_path))
        if not cands:
            print("[错误] 目录中未发现含 binder 链的复合物 PDB")
            sys.exit(1)
        print(f"[信息] 待评估复合物: {len(cands)} 个\n")
        rows = []
        for i, (sid, pdb_path) in enumerate(cands, 1):
            print(f"[{i}/{len(cands)}] 评估 {sid} ...")
            work_dir = os.path.join(tempfile.gettempdir(), f"schro_complex_{sid}")
            try:
                res = process_complex(pdb_path, args.target_chain, args.ligand_chain,
                                      work_dir, sid)
                e = res["energies"]
                dg = e.get("r_psp_MMGBSA_dG_Bind")
                dg_hbond = e.get("r_psp_MMGBSA_dG_Bind_Hbond")
                dg_vdw = e.get("r_psp_MMGBSA_dG_Bind_vdW")
                dg_solv = e.get("r_psp_MMGBSA_dG_Bind_Solv_GB")
                rows.append({
                    "id": sid,
                    "length": count_chain_residues(pdb_path, args.ligand_chain),
                    "mean_plddt": "",
                    "prime_dg": dg if dg is not None else "N/A",
                    "dg_hbond": dg_hbond if dg_hbond is not None else "",
                    "dg_vdw": dg_vdw if dg_vdw is not None else "",
                    "dg_solv_gb": dg_solv if dg_solv is not None else "",
                    "docked": 0,
                    "pass": "",
                    "note": "RFdiffusion复合物直算(无对接)",
                })
                print(f"  ΔG_bind = {dg if dg is not None else 'N/A'} kcal/mol")
            except Exception as e:
                print(f"  [错误] {sid}: {e}")
                rows.append({
                    "id": sid, "length": "", "mean_plddt": "", "prime_dg": "ERROR",
                    "dg_hbond": "", "dg_vdw": "", "dg_solv_gb": "", "docked": 0,
                    "pass": "", "note": f"计算失败: {e}",
                })
        if rows:
            fieldnames = ["id", "length", "mean_plddt", "prime_dg", "dg_hbond",
                          "dg_vdw", "dg_solv_gb", "docked", "pass", "note"]
            with open(args.out, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
                w.writeheader(); w.writerows(rows)
            print(f"\n[成功] 能量结果 -> {args.out}  ({len(rows)} 个复合物)")
        else:
            print("\n[警告] 无成功评估的复合物")
        return

    # 读取候选列表
    if args.candidates and os.path.exists(args.candidates):
        with open(args.candidates) as f:
            candidates = list(csv.DictReader(f))
    elif args.pdb_dir:
        candidates = []
        for pdb_path in glob.glob(os.path.join(args.pdb_dir, "*.pdb")):
            sid = os.path.splitext(os.path.basename(pdb_path))[0]
            candidates.append({"id": sid, "pass": "1"})
    else:
        print("[错误] 需提供 --candidates CSV 或 --pdb-dir")
        sys.exit(1)

    if args.pass_only:
        candidates = [c for c in candidates if c.get("pass") == "1"]
    if args.top_n:
        candidates = candidates[:args.top_n]

    print(f"[信息] 待评估候选: {len(candidates)} 个")
    print(f"[信息] 靶点: {args.target} (chain {args.target_chain})")
    if args.dock:
        print(f"[信息] 启用对接: prime_zbind (binder 链 {args.ligand_chain})")
    else:
        print(f"[信息] 未启用对接 —— 若 binder 为单独折叠，ΔG 可能无物理意义")

    rows = []
    for i, cand in enumerate(candidates, 1):
        sid = cand.get("id", f"cand_{i:03d}")
        pdb_path = os.path.join(args.pdb_dir or ".", f"{sid}.pdb")
        if not os.path.exists(pdb_path):
            print(f"  [跳过] {sid}: PDB 不存在 {pdb_path}")
            continue

        print(f"\n[{i}/{len(candidates)}] 评估 {sid} ...")
        work_dir = os.path.join(tempfile.gettempdir(), f"schrodinger_{sid}")

        try:
            res = process_candidate(
                pdb_path, args.target, args.target_chain, args.ligand_chain,
                work_dir, sid, args.dock,
            )
            energies = res["energies"]
            dg = energies.get("r_psp_MMGBSA_dG_Bind")
            dg_hbond = energies.get("r_psp_MMGBSA_dG_Bind_Hbond")
            dg_vdw = energies.get("r_psp_MMGBSA_dG_Bind_vdW")
            dg_solv = energies.get("r_psp_MMGBSA_dG_Bind_Solv_GB")

            note = "Schrödinger Prime MM-GBSA"
            if not res["docked"] and not args.dock:
                note += " | ⚠未对接:ΔG可能无物理意义"

            rows.append({
                "id": sid,
                "length": cand.get("length", ""),
                "mean_plddt": cand.get("mean_plddt", ""),
                "prime_dg": dg if dg is not None else "N/A",
                "dg_hbond": dg_hbond if dg_hbond is not None else "",
                "dg_vdw": dg_vdw if dg_vdw is not None else "",
                "dg_solv_gb": dg_solv if dg_solv is not None else "",
                "docked": 1 if res["docked"] else 0,
                "pass": cand.get("pass", ""),
                "note": note,
            })
            print(f"  ΔG_bind = {dg if dg is not None else 'N/A'} kcal/mol"
                  f"  (docked={res['docked']})")

        except Exception as e:
            print(f"  [错误] {sid}: {e}")
            rows.append({
                "id": sid,
                "length": cand.get("length", ""),
                "mean_plddt": cand.get("mean_plddt", ""),
                "prime_dg": "ERROR",
                "dg_hbond": "", "dg_vdw": "", "dg_solv_gb": "",
                "docked": 0, "pass": cand.get("pass", ""),
                "note": f"计算失败: {e}",
            })

    if rows:
        fieldnames = ["id", "length", "mean_plddt", "prime_dg", "dg_hbond",
                      "dg_vdw", "dg_solv_gb", "docked", "pass", "note"]
        with open(args.out, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
            w.writeheader(); w.writerows(rows)
        print(f"\n[成功] 能量结果 -> {args.out}  ({len(rows)} 个候选)")
    else:
        print("\n[警告] 无成功评估的候选")


if __name__ == "__main__":
    main()
