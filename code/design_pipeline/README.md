# code/design_pipeline

The **full de novo design pipeline**, provided for method transparency. Scripts are
grouped by stage. Several stages depend on **third-party licensed software**; those are
clearly marked and will not execute without the corresponding license. The
analysis / figure-generation layers (`../analysis`, `../figure_generation`) are fully
open and do not require these.

| Stage | Script(s) | External license required? |
|---|---|---|
| 0. Fetch target | `01_fetch_target_structure.py` | no (public PDB 1FLT) |
| 1. Hotspots | `02_define_hotspots.py` | no |
| 2. Backbone design | `02b_rfdiffusion_complex.py` | **yes — RFdiffusion commercial license** |
| 3. Sequence design | `02c_proteinmpnn_design.py` | no (ProteinMPNN open) |
| 4. Prescreen | `03_esmfold_prescreen.py` | no (ESMFold open) |
| 5. Re-rank | `04_boltz_rerank.py` | no (Boltz-1 open) |
| 5. Re-rank | `04_schrodinger_energy.py` | **yes — Schrödinger/Maestro commercial** |
| 5. Re-rank | `04_af3_rerank.py` | **yes — AlphaFold3 academic license** |
| 6. De-amyloid redesign | `05h_deA_redesign_1_2.py`, `05k_boltz_redesign_2.py` | no (Boltz-1) |
| 7. Relax + MM-GBSA | `05l_relax_mmgbsa_2deA.py`, `05m_interface_redesign_2.py` | **yes — Schrödinger/Maestro** (MM-GBSA step) |

> In the published work, **Boltz-1** was used as the AF3-class re-rank surrogate (no
> academic license needed) and **Schrödinger MM-GBSA** supplied the physical ΔG. The
> RFdiffusion and AF3 steps were run under their respective licenses.

## Inputs / outputs

- Target construct: 98-aa VEGF-A RBD (PDB 1FLT chain W, residues 1–98) — see
  `../../data/templates/`.
- Design outputs (sequences, PDBs, CRO constructs) live in the full project's
  `designs/` and `results/` directories; the 311 MB raw predictions are on Zenodo
  (see `../../results_manifest.md`).
