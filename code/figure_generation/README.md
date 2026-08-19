# code/figure_generation

Scripts that **regenerate every published figure** from the manuscript. All numeric
values used in the figures are hard-coded from the real pipeline outputs, so no 311 MB
download is required.

| Script | Produces | Depends on |
|---|---|---|
| `make_figures.py` | Figs. 1–8 (PNG + SVG) | NumPy, Matplotlib (self-contained data) |
| `make_graphical_abstract.py` | `Graphical_Abstract.png` | Matplotlib + the two PNGs below |
| `render_ga_outcome.py` | `GA_V1_contact.png`, `GA_V4_bridge.png` | **PyMOL** (open-source); reads `data/representative_structures/` |
| `make_fig7_pymol.py` | Fig. 7 panels (V1/V2/V4) | **PyMOL** (open-source); reads `data/representative_structures/` |
| `verify_figures.py` | guardrail: figures match manuscript text | reads `../../论文/` + regenerates Fig. 4 |

## Run

```bash
conda activate anti_vegf_repro        # or: pip install -r ../../requirements.txt
python make_figures.py                # Figs. 1–8

# Graphical Abstract (needs PyMOL)
/root/miniconda3/envs/pymol/bin/python render_ga_outcome.py
python make_graphical_abstract.py

# Fig. 7 (needs PyMOL)
/root/miniconda3/envs/pymol/bin/python make_fig7_pymol.py

# Guardrail
python verify_figures.py
```

## Notes

- `verify_figures.py` asserts the candidate→column mapping in Fig. 4 and that the Fig. 2 /
  Fig. 6 captions use the correct verbs ("improved" vs "maintained"; "both chains" vs
  "one chain"). Run it before committing any figure change.
- Paths are resolved **relative to this repository** (no absolute Windows paths), so the
  package runs wherever it is cloned.
