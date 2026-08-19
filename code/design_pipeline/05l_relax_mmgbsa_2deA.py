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
05l_relax_mmgbsa_2deA.py — #2 去A重设计 top 候选的 松弛 + MM-GBSA 双指标终评
对 05k 产出的"好模式" Boltz 复合物（seed 0/42，ipTM>0.5）做:
  PDBFixer 加氢 -> OpenMM Cα约束最小化(解界面互穿) -> Schrödinger prime_mmgbsa
与基线 #2 (vegf_len100_4, 原始 ΔG≈-37.8；de-A 目标: 低聚集且结合不劣化) 同口径。

仅用"好模式"PDB（排除 seed11 伪构象 ipTM≈0.10），与 T0.1_s35 3× 取好模式一致。
"""
import os, csv, glob, importlib.util

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, ".."))
SPEC = importlib.util.spec_from_file_location(
    "schro_energy", os.path.join(HERE, "04_schrodinger_energy.py"))
SCH = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SCH)

from pdbfixer import PDBFixer
from openmm import app, unit, CustomExternalForce, LangevinMiddleIntegrator

BOLTZ_PDB_DIR = os.path.join(PROJ, "results", "deA_redesign_1_2", "boltz_pdb")
WORK = os.path.join(PROJ, "results", "deA_redesign_1_2", "mmgbsa_work")
os.makedirs(WORK, exist_ok=True)
OUT_CSV = os.path.join(PROJ, "results", "deA_redesign_1_2", "mmgbsa_redesign_2.csv")
REF_DG_BASE = -37.8  # 原始 #2 (vegf_len100_4) 基线 ΔG 近似

# 好模式 PDB（ipTM>0.5，排除 seed11 伪构象）
GOOD = ["#2_T0.2_s83_s0_model_0.pdb",
        "#2_T0.2_s83_s42_model_0.pdb",
        "#2_T0.2_s68_s0_model_0.pdb",
        "#2_T0.2_s68_s42_model_0.pdb"]


def relax_one(inpdb, wd):
    os.makedirs(wd, exist_ok=True)
    heavy = os.path.join(wd, "heavy.pdb")
    with open(inpdb) as fin, open(heavy, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                if line[12:16].strip().startswith("H"):
                    continue
                fout.write(line)
            elif line.startswith(("TER", "END")):
                fout.write(line)
    fixer = PDBFixer(filename=heavy)
    fixer.findMissingResidues(); fixer.findMissingAtoms()
    fixer.addMissingAtoms(); fixer.addMissingHydrogens(7.0)
    fixed_pdb = os.path.join(wd, "fixed.pdb")
    app.PDBFile.writeFile(fixer.topology, fixer.positions, open(fixed_pdb, "w"))
    pdb = app.PDBFile(fixed_pdb)
    ff = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')
    system = ff.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic,
                              nonbondedCutoff=2.0 * unit.nanometer, constraints=None)
    k = 100.0 * unit.kilojoules_per_mole / (unit.nanometer ** 2)
    rest = CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    rest.addGlobalParameter('k', k)
    rest.addPerParticleParameter('x0'); rest.addPerParticleParameter('y0'); rest.addPerParticleParameter('z0')
    for i, at in enumerate(pdb.topology.atoms()):
        if at.name == 'CA':
            p = pdb.positions[i]
            rest.addParticle(i, [p[0].value_in_unit(unit.nanometer),
                                 p[1].value_in_unit(unit.nanometer),
                                 p[2].value_in_unit(unit.nanometer)])
    system.addForce(rest)
    integ = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.004 * unit.picosecond)
    sim = app.Simulation(pdb.topology, system, integ)
    sim.context.setPositions(pdb.positions)
    sim.minimizeEnergy(maxIterations=8000, tolerance=5.0 * unit.kilojoule_per_mole / unit.nanometer)
    minpdb = os.path.join(wd, "minimized_complex.pdb")
    app.PDBFile.writeFile(pdb.topology, sim.context.getState(getPositions=True).getPositions(), open(minpdb, 'w'))
    return minpdb


def main():
    if not SCH.check_schrodinger():
        sys.exit(1)
    rows = []
    for pdb_name in GOOD:
        inpdb = os.path.join(BOLTZ_PDB_DIR, pdb_name)
        if not os.path.exists(inpdb):
            print(f"[跳过] 缺 {pdb_name}")
            continue
        job = pdb_name.replace("_model_0.pdb", "")
        print(f"\n=== {job} ===", flush=True)
        wd = os.path.join(WORK, job)
        minpdb = relax_one(inpdb, wd)
        try:
            res = SCH.process_complex(minpdb, "A", "B", os.path.join(wd, "energy"), job)
            e = res["energies"]
            dg = e.get("r_psp_MMGBSA_dG_Bind")
            print(f"  ΔG_bind = {dg} kcal/mol (原始#2基线≈{REF_DG_BASE})", flush=True)
            rows.append({"id": job, "deltaG_kcal": dg,
                         "note": "PDBFixer+OpenMM Cα约束最小化+MM-GBSA(好模式)"})
        except Exception as ex:
            print(f"  [MM-GBSA失败] {ex}", flush=True)
            rows.append({"id": job, "deltaG_kcal": "ERROR", "note": f"失败:{ex}"})
    if rows:
        with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["id", "deltaG_kcal", "note"])
            w.writeheader(); w.writerows(rows)
        # 汇总每候选均值
        from collections import defaultdict
        agg = defaultdict(list)
        for r in rows:
            if isinstance(r["deltaG_kcal"], (int, float)):
                agg[r["id"].rsplit("_s", 1)[0]].append(r["deltaG_kcal"])
        print("\n=== 每候选 MM-GBSA ΔG 均值 ===")
        for cand, vals in agg.items():
            m = round(sum(vals) / len(vals), 1)
            print(f"  {cand}: {vals} -> 均值 {m} kcal/mol  (vs 原始#2 {REF_DG_BASE})")
        print(f"\n[成功] -> {OUT_CSV}")


if __name__ == "__main__":
    main()
