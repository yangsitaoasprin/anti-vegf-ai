# code/analysis

Scripts that compute the **reported scores and tables** from the pipeline outputs. These
run on open-source tooling (Boltz-1, OpenFold3, NumPy, etc.) and read the raw prediction
directories (`results/...` in the full project) or the bundled tables in
`../../data/manuscript_tables/`.

| Script | Role in the paper |
|---|---|
| `04_boltz_rerank.py` | Boltz-1 ipTM / interface-pAE re-ranking (AF3-class surrogate) |
| `04_schrodinger_energy.py` | Schrödinger MM-GBSA ΔG scoring (physical re-rank) |
| `04b_boltz_comparison.py` | Boltz vs AF3-style score comparison |
| `05_rank_and_report.py`, `05_rank_candidates.py` | composite ranking → final candidate set |
| `05f_mmgbsa_redesign.py` | de-amyloid redesign MM-GBSA re-evaluation (core result) |
| `05x_boltz_fullVEGF.py` | **full-length VEGF-A165 re-check** (validates truncation bias) |
| `06_of3_crossval.py` | OpenFold3 cross-validation → Fig. 4 |
| `05b_waltz_aggregation.py` | WALTZ aggregation-risk filter |
| `05o_mmgbsa_if.py`, `05r_energy_breakdown.py`, `05w_parse_conservative.py` | MM-GBSA interface / energy decomposition / conservative-mutation parsing |
| `analyze_binding_mode.py`, `analyze_binding_mode_4.py`, `analyze_binding_mode_5.py` | per-candidate binding-mode analysis |

## Reproducibility note

- The per-seed MM-GBSA decomposition that backs **Table S2** is bundled as
  `../../data/manuscript_tables/tableS2_perseed_mmgbsa.csv`.
- The composite final ranking (Table 2 / candidate set) is bundled as
  `../../data/manuscript_tables/candidates_6_final.csv`.
- Re-running these scripts end-to-end requires the raw 311 MB model outputs (Zenodo DOI
  in `../../results_manifest.md`); the bundled CSVs are sufficient to inspect the numbers.
