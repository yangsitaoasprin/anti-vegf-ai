# Reproducible code & data — *Developability-gated virtual screening of de novo anti-VEGF-A miniprotein binders*

This repository provides the code and minimal supporting data needed to reproduce the
figures, tables, and physical-ranking results reported in the accompanying CSBJ
manuscript (Manuscript_v7). The work designs and ranks de novo miniprotein binders
against VEGF-A (UniProt P15692-4) using a **developability-gated, multi-model virtual
screening pipeline**:

1. **Design layer** — RFdiffusion (complex design) → ProteinMPNN (sequence) →
   ESMFold (prescreen) → Boltz-1 (re-rank, AF3-class local surrogate) →
   Schrödinger MM-GBSA ΔG re-rank.
2. **Adjudication layer** — an MM-GBSA van-der-Waals (vdW) gate that flags
   steric-clash poses, plus a template-guided OpenFold3 (OF3) cross-validation.
3. **MD cross-check layer (§2.10)** — 5-ns restrained implicit-solvent molecular
   dynamics (OpenMM 8.4) per final candidate, with per-frame single-point MM-GBSA
   re-scoring, confirming the single-conformation ranking is not a frozen-pose artifact.

## Core contributions (aligned with Manuscript_v7)

| # | Contribution | Validation |
|---|--------------|------------|
| ① | **MM-GBSA vdW gate** — binary rule `vdW < 0 AND ΔG < 0 = clean binding` rejects ipTM-favoured but geometrically clashing poses | single-structure MM-GBSA (vdW term) vs ipTM; V2 original pose vdW = +4527 kcal/mol at ipTM 0.765 |
| ② | **1FLT template-guided OF3** — a single high-resolution crystal template raises OF3 ipTM ~3.5× | seq-only vs 1FLT-tmpl vs ColabFold MSA (0.07–0.35, unstable) |
| ③ | **§2.10 GPU MD + MM-GBSA cross-check** — 5-ns restrained MD + per-frame MM-GBSA | OpenMM 8.4 / RTX 5090; per-frame ΔG ± SD |
| ④ | **Developability-driven redesign** — WALTZ + official APR-Score as two independent aggregation opinions; iterative de-alanine / de-amyloid / conservative tuning | pre/post redesign vdW-clean rate, APR window count, WALTZ pass rate |
| ⑤ | **V4 dimer-bridge hypothesis** — candidate-specific conformation contacting both VEGF-A chains | chain-pair ipTM ≈ 0.66–0.68 |

## Final candidates (delivery set)

| Cand. | ID | Spec | ipTM | ΔG (kcal/mol) |
|-------|----|------|------|---------------|
| V1 | `#1_T0.3_s132` | 80 aa, de-alanine | 0.937 | −130 |
| V2 | `#2c_T0.1_s18` | 100 aa, dual-receptor, A85E+K95A conservative tuning | 0.63 | −53 |
| V4 | `T0.1_s35` | 80 aa, de-amyloid 16 muts, WALTZ clean | 0.70 | −51 |

§2.10 per-frame MM-GBSA (kcal/mol): **V1 −74.1 ± 10.6, V2 −89.1 ± 11.6, V4 −53.3 ± 5.5**,
consistent with the single-conformation ranking (V1 ≥ V4 > V2).

## Author & contact

| | |
|---|---|
| **Author** | Yang Sitao (杨四涛) |
| **Affiliation** | College of Pharmacy, Dali University; The Third Affiliated Hospital of Dali University |
| **Email** | yangsitaoasprin@swpu.edu.cn |
| **Code license** | MIT (see `LICENSE`) |
| **Data license** | CC-BY-4.0 (representative structures + tables in `data/`); raw 311 MB model outputs under the Zenodo DOI below |

> Full author list, CRediT roles, competing interests and funding are declared in the
> manuscript. See `CODE_AUTHORS.md` for the contributor/contact block used in this repository.

## Repository layout

```
repro/
├── README.md                      # this file
├── LICENSE                        # MIT
├── CODE_AUTHORS.md                # author + contact + CRediT block
├── environment.yml                # conda env spec (analysis / figure generation)
├── requirements.txt               # pip-equivalent for the above
├── results_manifest.md            # inventory of the 311 MB raw model outputs (Zenodo)
├── UPLOAD_GUIDE.md                # how to push to GitHub + deposit raw data on Zenodo
├── code/
│   ├── figure_generation/         # scripts that regenerate Figs. 1–8 + Graphical Abstract
│   ├── analysis/                  # scripts that compute the reported scores / tables
│   └── design_pipeline/           # design + re-rank scripts (some need external licenses)
└── data/
    ├── representative_structures/ # 5 complex structures used for Fig. 7 & Graphical Abstract
    ├── templates/                 # 1FLT crystal template (PDB + mmCIF)
    └── manuscript_tables/         # candidates_6_final.csv + per-seed MM-GBSA (Table S2)
```

> **Note on the §2.10 MD cross-check:** the MD production/evaluation source
> (`06a_implicit_md.py`, `06b_trajectory_mmgbsa.py`) lives in the main project's
> `脚本/` directory, not bundled here. To make this repository fully self-contained for
> the MD layer, copy those two scripts (plus `openmm`/`pdbfixer`) into `code/analysis/`
> — see "Planned additions" below. The figure/table layers below are fully runnable as-is.

## Quick start (regenerate the published figures)

The figure-generation scripts are **self-contained**: all numeric values used in the
figures are hard-coded from the real pipeline outputs, so no 311 MB download is needed
to reproduce the figures.

```bash
# 1) (optional) create the conda env
conda env create -f environment.yml
conda activate anti_vegf_repro

# 2) regenerate Figs. 1–8 (writes PNG+SVG next to make_figures.py)
cd code/figure_generation
python make_figures.py

# 3) regenerate the Graphical Abstract panels (needs PyMOL; see below)
/root/miniconda3/envs/pymol/bin/python render_ga_outcome.py   # writes GA_V1_contact.png, GA_V4_bridge.png
python make_graphical_abstract.py                              # assembles Graphical_Abstract.png

# 4) guardrail: confirm figures still match the manuscript text
python verify_figures.py
```

`make_fig7_pymol.py` (Fig. 7) and `render_ga_outcome.py` (Graphical Abstract panels)
require **PyMOL** (open-source build, `conda install -c conda-forge pymol-open-source`).
They read only the small representative structures in `data/representative_structures/`.

> **Fig. 3 interface-size stats:** `make_fig8_interface_stats.py` is bundled here and
> recomputes the interface-residue counts directly from the three representative complexes
> in `data/representative_structures/` (5-angstrom heavy-atom cutoff, the rule stated in the
> Fig. 3 caption). It writes PNG + SVG + `Fig8_interface_stats_numbers.json`. The counts it
> returns -- V1 18/25 (31.2%), V2 27/29 (29.0%), V4 19/17 (21.2%) -- are the ones printed in
> the manuscript, so Fig. 3 reproduces from a clone with no raw-data download.
>
> **Supplementary Fig. S1** (APR-Score recovery) is regenerated by
> `make_figS1_apr_recovery.py`, likewise self-contained.

## code/analysis

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

> The §2.10 MD cross-check (`06a_implicit_md.py`, `06b_trajectory_mmgbsa.py`) is **not**
> bundled in this repository yet — see "Planned additions". Its results (per-frame ΔG above)
> are reported from the main project's `脚本/` runs.

## code/design_pipeline

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

## code/figure_generation

Scripts that **regenerate every published figure** from the manuscript. All numeric
values used in the figures are hard-coded from the real pipeline outputs, so no 311 MB
download is required.

| Script | Produces | Depends on |
|---|---|---|
| `make_figures.py` | Figs. 1–8 (PNG + SVG) | NumPy, Matplotlib (self-contained data) |
| `make_fig8_interface_stats.py` | Fig. 3 interface counts (PNG + SVG + JSON) | NumPy, Matplotlib; reads `data/representative_structures/` |
| `make_figS1_apr_recovery.py` | Supplementary Fig. S1 (PNG + SVG) | NumPy, Matplotlib (self-contained data) |
| `make_graphical_abstract.py` | `Graphical_Abstract.png` (+ PDF/TIFF) | Matplotlib only -- the PyMOL panels below are rendered for the record but are **not** embedded in the Graphical Abstract |
| `render_ga_outcome.py` | `GA_V1_contact.png`, `GA_V4_bridge.png` | **PyMOL** (open-source); reads `data/representative_structures/` |
| `make_fig7_pymol.py` | Fig. 7 panels (V1/V2/V4) | **PyMOL** (open-source); reads `data/representative_structures/` |
| `figstyle.py` | *library*: enforces Elsevier's 7 pt printed-text minimum | Matplotlib; reads `FIG_FILES`/`FIG_WIDTH` from `build_docx.py` when that file is present, otherwise its built-in width table |
| `verify_figures.py` | guardrail: figures match manuscript text | **authoring tree only** -- needs the `论文/` folder, which is not shipped here |

- `verify_figures.py` asserts the candidate→column mapping in Fig. 4 and that the Fig. 2 /
  Fig. 6 captions use the correct verbs ("improved" vs "maintained"; "both chains" vs
  "one chain"). Run it before committing any figure change.
- Paths are resolved **relative to this repository** (no absolute Windows paths), so the
  package runs wherever it is cloned.

## Notes on scope & reproducibility

- **Minimal supporting data.** To keep this repository light, only the representative
  complex structures and the final/per-seed tables are bundled. The full 311 MB of raw
  Boltz-1 (`boltz_pdb/`) and OpenFold3 (`of3_out*/`) predictions — from which every number
  in the paper is derived — are deposited on **Zenodo** (DOI placeholder in `results_manifest.md`). Reviewers can therefore reproduce any individual pose.
- **Aggregation-scoring record.** The APR-Score / CCS recovery chain behind the
  aggregation screening and Supplementary Fig. S1 — scripts, derived tables and
  internal reports — has its own **Zenodo** record:
  `https://doi.org/10.5281/zenodo.ZZZZZZZZ`.
- **Design-generation scripts.** `code/design_pipeline/` contains the full design pipeline,
  including steps that require third-party licenses (RFdiffusion commercial license,
  Schrödinger/Maestro commercial license, AlphaFold3 academic license). Those scripts are
  provided for **method transparency**; they are clearly marked and will not run without
  the corresponding licensed software. The analysis and figure-generation layers run
  entirely on open tools (Boltz-1, OpenFold3, PyMOL open-source, NumPy, Matplotlib).
- **Target construct.** All design work targets the 98-aa VEGF-A RBD (PDB 1FLT chain W,
  residues 1–98). Wet-lab validation must use full-length VEGF-A165 (UniProt P15692-4,
  includes the heparin-binding domain 111–165).

## Planned additions (to make this repository self-contained)

- Bundle `06a_implicit_md.py` + `06b_trajectory_mmgbsa.py` from the main project `脚本/`
  into `code/analysis/`, so the §2.10 MD cross-check is reproducible from this repository.
*Done:* `make_fig8_interface_stats.py`, `make_figS1_apr_recovery.py` and `figstyle.py` are
now bundled in `code/figure_generation/`, so the entire figure layer reproduces from a
clone. Verified by copying the bundle into a directory with no `论文/` folder and running
every generator there (all exit 0, nothing written outside the clone).

## Citation

If you use this code or data, please cite the manuscript and the Zenodo record (DOI to be added after deposit).
