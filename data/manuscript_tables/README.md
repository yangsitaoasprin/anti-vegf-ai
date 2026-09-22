# Supplementary tables - machine-readable backing files

These CSV files carry the same values as the supplementary tables of the manuscript
*Developability-gated virtual screening of de novo anti-VEGF-A miniprotein binders*.
S1, S3, S4 and S5 were extracted mechanically from the manuscript source; no value was
recomputed, rounded or edited by hand. Where the manuscript prints an en-dash range
(e.g. `49-54`) or a Unicode minus (U+2212), the CSV keeps the manuscript character
verbatim.

| File | Manuscript table | Content |
|---|---|---|
| `tableS1_waltz_validation.csv` | Table S1 | WALTZ validation on the APR balanced test sample: TP/TN, precision, sensitivity |
| `tableS2_perseed_mmgbsa.csv` | Table S2 | Full-length VEGF-A165 per-seed MM-GBSA decomposition |
| `tableS3_seed_sensitivity.csv` | Table S3 | Full-length Boltz-1 ipTM seed sensitivity, 3-seed vs 5-seed panels |
| `tableS4_patch_summary.csv` | Table S4 | BioLuminate Patch Analyzer per-construct summary (free-state monomer) |
| `tableS5_patch_windows.csv` | Table S5 | All aggregation-risk windows individually, with free/bound accessibility |
| `candidates_6_final_pre_redesign.csv` | Table S6 | Original Top-6 designs selected by the composite score (pre-redesign) |
| `candidates_6_final.csv` | Table 1 | Final candidate table for the six redesigns (post-redesign) |

Raw inputs behind these tables are in `../per_seed/` (Tables S2/S3) and
`../patch_analyzer/` (Tables S4/S5).
