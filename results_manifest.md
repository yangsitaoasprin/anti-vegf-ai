# Raw model-output manifest (Zenodo deposit)

The published figures and tables are fully reproducible from the **bundled** data
(`data/`) plus the hard-coded values in `code/figure_generation/`. For reviewers who want
to re-derive **any individual pose** or re-run the scoring, the complete raw predictions
are deposited on **Zenodo** (CC-BY-4.0).

> **Zenodo DOI (placeholder — fill after upload):** `10.5281/zenodo.XXXXXXX`
> Upload steps: see `UPLOAD_GUIDE.md`.
> **Related deposit — APR-Score / CCS recovery record:** `https://doi.org/10.5281/zenodo.ZZZZZZZZ`
> (the scripts, derived tables and internal reports behind the aggregation
> screening and Supplementary Fig. S1).

## What is in the Zenodo archive (≈ 311 MB)

Root: `anti_VEGF_AI_project/results/` (only the model-output directories are archived;
intermediate logs and deprecated files are excluded).

| Path | Size | Content |
|---|---|---|
| `deA_redesign_1_2/of3_out/` | ~128 MB | OpenFold3 predictions, 98-aa RBD construct (mmCIF + JSON scores) |
| `deA_redesign_1_2/of3_out_std165/` | ~78 MB | OpenFold3 predictions, full-length VEGF-A165 (incl. template & dimer modes) |
| `deA_redesign_1_2/boltz_pdb/` | ~38 MB | Boltz-1 predicted complexes (PDB) for all candidates/seeds |
| `deA_redesign_1_2/boltz_pdb_fullVEGF/` | ~5 MB | Boltz-1 full-length VEGF-A165 complexes |
| `deA_redesign_1_2/mmgbsa_work_gen/` | ~31 MB | Schrödinger MM-GBSA working files (per pose) |
| `deA_redesign_1_2/mmgbsa_work*/` | ~52 MB | additional MM-GBSA working dirs |
| `results/mmgbsa_redesign/` | ~8.6 MB | de-amyloid redesign MM-GBSA inputs/PD bs |
| `deA_redesign_1_2/of3_pdb*/`, `esmfold_pdb*/`, `clash_sweep/` | ~4 MB | supporting structures & clash sweeps |
| `deA_redesign_1_2/*.csv`, `*.md`, `fullVEGF_*.csv` | < 1 MB | per-seed score tables, full-length re-check reports |

## Provenance of the key numbers

| Manuscript element | Source in archive |
|---|---|
| Table 2 (composite ranking) | `results/candidates_6_final.csv` (also bundled in `data/manuscript_tables/`) |
| Table S2 (per-seed MM-GBSA) | `deA_redesign_1_2/fullVEGF_mmgbsa_3seed.csv` (+ bundled `tableS2_perseed_mmgbsa.csv`) |
| Fig. 4 cross-validation | `deA_redesign_1_2/of3_iptm.csv`, `of3_tmpl_iptm.csv`, `of3_mmgbsa.csv`, `of3_tmpl_mmgbsa.csv` |
| Fig. 7 / Graphical Abstract | representative structures bundled in `data/representative_structures/` (derived from `boltz_pdb/` + `of3_out_std165/`) |

## Not archived (by design)

- RFdiffusion / ProteinMPNN intermediate trajectories (regenerable from seeds; very large).
- Deprecated intermediates (`results/DEPRECATED_INTERMEDIATES.md`, `*.DEPRECATED*` files).
- Local licensed-software caches (Schrödinger/PyRosetta) — not redistributable.
