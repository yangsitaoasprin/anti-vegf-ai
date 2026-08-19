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
05f_mmgbsa_redesign.py
对 Plan A ProteinMPNN 重设计的 top 候选 (Boltz 复合物 PDB) 批量跑 Schrödinger
Prime MM-GBSA，与原始 #4 (vegf_len80_3, ΔG≈-37.8 kcal/mol) 对比。

输入: results/mmgbsa_redesign/pdbs/*.pdb  (Boltz 预测复合物, 链 A=VEGF 靶点, B=微蛋白)
输出: results/mmgbsa_redesign/mmgbsa_results.csv
"""
import csv
import glob
import os
import sys
import io
import importlib.util
import argparse

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(SCRIPT_DIR)

# 数字开头的模块名无法直接 import，用 importlib 按路径加载
_spec = importlib.util.spec_from_file_location(
    "schrodinger_energy", os.path.join(SCRIPT_DIR, "04_schrodinger_energy.py")
)
se = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(se)

PDB_DIR = os.path.join(PROJ, "results", "mmgbsa_redesign", "pdbs")
OUT_CSV = os.path.join(PROJ, "results", "mmgbsa_redesign", "mmgbsa_results.csv")
WORK_ROOT = os.path.join(PROJ, "results", "mmgbsa_redesign", "work")
TARGET_CHAIN = "A"
LIGAND_CHAIN = "B"

REF_DG_WT = -37.8  # 原始 #4 的 MM-GBSA ΔG 基线 (kcal/mol)


def main():
    ap = argparse.ArgumentParser(description="MM-GBSA 批量评估 Boltz 复合物")
    ap.add_argument("--only", default=None, help="只处理文件名包含此子串的 PDB（如 T0.1_s35）")
    args = ap.parse_args()
    if not se.check_schrodinger():
        sys.exit(1)
    pdbs = sorted(glob.glob(os.path.join(PDB_DIR, "*.pdb")))
    if args.only:
        pdbs = [p for p in pdbs if args.only in os.path.basename(p)]
    if not pdbs:
        print(f"[错误] {PDB_DIR} 下无 PDB")
        sys.exit(1)
    os.makedirs(WORK_ROOT, exist_ok=True)
    print(f"[信息] 待评估复合物: {len(pdbs)} 个 (filter={args.only})")
    rows = []
    for i, pdb in enumerate(pdbs, 1):
        sid = os.path.splitext(os.path.basename(pdb))[0]
        job = sid.replace("_model_0", "")
        print(f"\n[{i}/{len(pdbs)}] {job}  ...")
        work_dir = os.path.join(WORK_ROOT, job)
        try:
            res = se.process_complex(pdb, TARGET_CHAIN, LIGAND_CHAIN, work_dir, job)
            e = res["energies"]
            dg = e.get("r_psp_MMGBSA_dG_Bind")
            rows.append({
                "id": job,
                "prime_dg": dg if dg is not None else "N/A",
                "dg_hbond": e.get("r_psp_MMGBSA_dG_Bind_Hbond", ""),
                "dg_vdw": e.get("r_psp_MMGBSA_dG_Bind_vdW", ""),
                "dg_solv_gb": e.get("r_psp_MMGBSA_dG_Bind_Solv_GB", ""),
                "dg_coulomb": e.get("r_psp_MMGBSA_dG_Bind_Coulomb", ""),
                "dg_covalent": e.get("r_psp_MMGBSA_dG_Bind_Covalent", ""),
                "dg_lipo": e.get("r_psp_MMGBSA_dG_Bind_Lipo", ""),
                "dg_packing": e.get("r_psp_MMGBSA_dG_Bind_Packing", ""),
                "dg_self_cont": e.get("r_psp_MMGBSA_dG_Bind_SelfCont", ""),
                "note": "Boltz复合物直算(无对接)",
            })
            print(f"  ΔG_bind = {dg if dg is not None else 'N/A'} kcal/mol  "
                  f"(原始#4 基线: {REF_DG_WT})")
        except Exception as ex:
            print(f"  [错误] {job}: {ex}")
            rows.append({
                "id": job, "prime_dg": "ERROR", "dg_hbond": "", "dg_vdw": "",
                "dg_solv_gb": "", "dg_coulomb": "", "dg_covalent": "",
                "dg_lipo": "", "dg_packing": "", "dg_self_cont": "",
                "note": f"计算失败: {ex}",
            })
    if rows:
        fieldnames = ["id", "prime_dg", "dg_hbond", "dg_vdw", "dg_solv_gb",
                      "dg_coulomb", "dg_covalent", "dg_lipo", "dg_packing",
                      "dg_self_cont", "note"]
        with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        print(f"\n[成功] 能量结果 -> {OUT_CSV}  ({len(rows)} 个复合物)")
    else:
        print("\n[警告] 无成功评估的复合物")


if __name__ == "__main__":
    main()
