# Upload guide — GitHub (code) + Zenodo (raw data)

This repository is **local only** so far. Two public deposits are required for CSBJ
compliance (Elsevier does **not** accept "available on request").

> You (the author) must perform the actual push / upload — they need your credentials.
> This guide tells you exactly what to do. After each deposit, paste the resulting DOI
> back into `README.md`, `results_manifest.md`, `CODE_AUTHORS.md`, and the manuscript's
> Data/Code availability statement.

---

## 1. GitHub — the code & minimal data in this `repro/` folder

```bash
cd /path/to/anti_VEGF_AI_project/repro
git init
git add -A
git commit -m "Reproducible code & minimal data for anti-VEGF-A miniprotein screen"

# Create the repo on GitHub (e.g. via https://github.com/new or gh CLI), then:
git branch -M main
git remote add origin https://github.com/yangsitaoasprin/<repo-name>.git
git push -u origin main
```

- The `repro/` folder is small (~1.4 MB of bundled data + scripts) — well within GitHub
  limits; **do not** `git add` the 311 MB `results/` raw outputs here.
- Get the **GitHub DOI** via the GitHub "Create release" → "Publish" flow (issues a
  Zenodo/GitHub DOI), or simply use the repo URL + a tagged release.

---

## 2. Zenodo — the 311 MB raw model outputs

1. Log in to https://zenodo.org (preferably via GitHub single sign-on).
2. **Upload** → new deposit.
3. Drag in `anti_VEGF_AI_project/results/` **excluding** logs/deprecated files (see
   `results_manifest.md` for the exact directory list). Or compress first:
   ```bash
   cd /path/to/anti_VEGF_AI_project
   tar czf anti_vegf_raw_results.tar.gz \
     results/deA_redesign_1_2 results/mmgbsa_redesign results/candidates_6_final.csv
   ```
4. Metadata:
   - **Title**: *Raw model outputs for "Developability-gated, multi-model virtual screening of de novo anti-VEGF-A miniprotein binders"*
   - **Authors**: Yang Sitao (affiliation: Dali University)
   - **License**: Creative Commons Attribution 4.0 (CC-BY-4.0)
   - **Resource type**: Dataset
   - **Related identifiers**: link the GitHub repo (relation: "is supplemented by").
5. Publish → copy the **DOI** (e.g. `10.5281/zenodo.XXXXXXX`).

---

## 3. Backfill the DOIs

Replace the placeholders:

| File | Placeholder to replace |
|---|---|
| `README.md` | Zenodo DOI line |
| `results_manifest.md` | `10.5281/zenodo.XXXXXXX` |
| `CODE_AUTHORS.md` | `<zenodo-DOI>`, `<repo-URL>` |
| `论文/Manuscript_v7.md` (Data/Code availability) | GitHub URL + Zenodo DOI |
| `论文/Manuscript_ZH.md` (数据/代码可用性) | 同上（中文） |

---

## 4. CSBJ Data/Code availability — required wording (already drafted in the manuscript)

> *Data availability.* Representative complex structures, the template (PDB 1FLT), the
> final candidate table and the per-seed MM-GBSA decomposition are available in the
> GitHub repository at `<repo-URL>` (DOI: `<github-DOI>`). The complete raw Boltz-1 and
> OpenFold3 predictions (~311 MB) are deposited on Zenodo at
> `https://doi.org/<zenodo-DOI>` under CC-BY-4.0. VEGF-A structural data derive from PDB
> entry 1FLT. Code to regenerate all figures and tables is in the same GitHub repository
> (MIT license).
