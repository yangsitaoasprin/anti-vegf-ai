# Reproducible code & data — *Developability-gated, multi-model virtual screening of de novo anti-VEGF-A miniprotein binders*

This repository provides the code and minimal supporting data needed to reproduce the
figures, tables, and physical-ranking results reported in the accompanying CSBJ
manuscript. The work designs and ranks de novo miniprotein binders against VEGF-A
(UniProt P15692-4) using a **two-layer virtual-screening pipeline**:

1. **Design layer** — RFdiffusion (complex design) → ProteinMPNN (sequence) →
   ESMFold (prescreen) → Boltz-1 (re-rank, AF3-class local surrogate) →
   Schrödinger MM-GBSA ΔG re-rank.
2. **Adjudication layer** — an MM-GBSA van-der-Waals (vdW) gate that flags
   steric-clash poses, plus a template-guided OpenFold3 (OF3) cross-validation.

## Author & contact

| | |
|---|---|
| **Author** | Yang Sitao (杨四涛) |
| **Affiliation** | School of Pharmaceutical Sciences, Dali University; The Third Affiliated Hospital of Dali University |
| **Email** | yangsitaoasprin@126.com |
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

## Notes on scope & reproducibility

- **Minimal supporting data.** To keep this repository light, only the representative
  complex structures and the final/per-seed tables are bundled. The full 311 MB of raw
  Boltz-1 (`boltz_pdb/`) and OpenFold3 (`of3_out*/`) predictions — from which every number
  in the paper is derived — are deposited on **Zenodo** (DOI placeholder in `results_manifest.md`). Reviewers can therefore reproduce any individual pose.
- **Design-generation scripts.** `code/design_pipeline/` contains the full design pipeline,
  including steps that require third-party licenses (RFdiffusion commercial license,
  Schrödinger/Maestro commercial license, AlphaFold3 academic license). Those scripts are
  provided for **method transparency**; they are clearly marked and will not run without
  the corresponding licensed software. The analysis and figure-generation layers run
  entirely on open tools (Boltz-1, OpenFold3, PyMOL open-source, NumPy, Matplotlib).
- **Target construct.** All design work targets the 98-aa VEGF-A RBD (PDB 1FLT chain W,
  residues 1–98). Wet-lab validation must use full-length VEGF-A165 (UniProt P15692-4,
  includes the heparin-binding domain 111–165).

## Citation

If you use this code or data, please cite the manuscript and the Zenodo record (DOI to be added after deposit).
