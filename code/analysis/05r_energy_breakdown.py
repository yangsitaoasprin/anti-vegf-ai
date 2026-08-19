# =============================================================================
# Author : Yang Sitao (杨四涛)
# Affiliation : School of Pharmaceutical Sciences, Dali University;
#               The Third Affiliated Hospital of Dali University
# Contact : yangsitaoasprin@126.com
# License : MIT (see ../../LICENSE)
# =============================================================================
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""拆解 MM-GBSA 能量项：对比 #2if s12 与 s83 的 E_complex/E_receptor/E_ligand，
定位 +9500 爆炸项（是复合物内部、还是受体/配体内部）。"""
import os, sys, importlib.util
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, ".."))
SPEC = importlib.util.spec_from_file_location("schro_energy", os.path.join(HERE, "04_schrodinger_energy.py"))
SCH = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(SCH)
from pdbfixer import PDBFixer
from openmm import app, unit, CustomExternalForce, LangevinMiddleIntegrator

BOLTZ = os.path.join(PROJ, "results", "deA_redesign_1_2", "boltz_pdb")
WORK = os.path.join(PROJ, "results", "deA_redesign_1_2", "mmgbsa_diag")
os.makedirs(WORK, exist_ok=True)

def relax(inpdb, outpdb):
    fixer = PDBFixer(filename=inpdb)
    fixer.findMissingResidues(); fixer.findMissingAtoms()
    fixer.addMissingAtoms(); fixer.addMissingHydrogens(7.0)
    fixed = os.path.join(WORK, "fixed.pdb")
    app.PDBFile.writeFile(fixer.topology, fixer.positions, open(fixed, "w"))
    pdb = app.PDBFile(fixed)
    ff = app.ForceField('amber99sbildn.xml', 'amber99_obc.xml')
    system = ff.createSystem(pdb.topology, nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=2.0*unit.nanometer, constraints=None)
    k = 100.0*unit.kilojoules_per_mole/(unit.nanometer**2)
    rest = CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    rest.addGlobalParameter('k', k)
    rest.addPerParticleParameter('x0'); rest.addPerParticleParameter('y0'); rest.addPerParticleParameter('z0')
    for i, at in enumerate(pdb.topology.atoms()):
        if at.name == 'CA':
            p = pdb.positions[i]
            rest.addParticle(i, [p[0].value_in_unit(unit.nanometer), p[1].value_in_unit(unit.nanometer), p[2].value_in_unit(unit.nanometer)])
    system.addForce(rest)
    integ = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.004*unit.picosecond)
    sim = app.Simulation(pdb.topology, system, integ)
    sim.context.setPositions(pdb.positions)
    sim.minimizeEnergy(maxIterations=8000, tolerance=5.0*unit.kilojoule_per_mole/unit.nanometer)
    app.PDBFile.writeFile(pdb.topology, sim.context.getState(getPositions=True).getPositions(), open(outpdb, "w"))

for label, pdb in [("#2if_s12_s0", os.path.join(BOLTZ, "#2if_T0.1_s12_s0_model_0.pdb")),
                   ("deA_s83_s0", os.path.join(BOLTZ, "#2_T0.2_s83_s0_model_0.pdb"))]:
    out = os.path.join(WORK, label + "_min.pdb")
    relax(pdb, out)
    res = SCH.process_complex(out, "A", "B", os.path.join(WORK, "en_" + label), label)
    e = res["energies"]
    print(f"\n=== {label} ===")
    for k, v in e.items():
        print(f"  {k}: {v}")
