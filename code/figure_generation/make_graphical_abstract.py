# =============================================================================
# Author : Yang Sitao (杨四涛)
# Affiliation : School of Pharmaceutical Sciences, Dali University;
#               The Third Affiliated Hospital of Dali University
# Contact : yangsitaoasprin@swpu.edu.cn
# License : MIT (see ../../LICENSE)
# Purpose: Assembles the Graphical Abstract (matplotlib) using the two PyMOL
#          panels produced by render_ga_outcome.py.
# =============================================================================
"""Graphical Abstract for CSBJ submission -- self-contained vector schematic.

Single coordinate system: main ax 0..10 x 0..13. All shapes drawn in data coords.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle
from matplotlib.image import imread

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(_HERE, "Graphical_Abstract.png")
C_V1 = "#2a6f97"; C_V2 = "#c1121f"; C_V4 = "#2a9d8f"
C_ACC = "#e9c46a"; C_GREY = "#6c757d"; C_NAVY = "#0d3a5c"
plt.rcParams["font.family"] = "DejaVu Sans"

fig = plt.figure(figsize=(10, 13), dpi=150)
ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
ax.set_xlim(0, 10); ax.set_ylim(0, 13)

def box(x, y, w, h, fc, ec="#333", lw=1.2, r=0.06, style="round,pad=0.02,rounding_size=0.05"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=style, fc=fc, ec=ec, lw=lw))

def T(x, y, s, sz=11, col="#1a1a1a", ha="left", va="bottom", wt="normal", st="normal", bbox=None):
    ax.text(x, y, s, fontsize=sz, color=col, ha=ha, va=va, weight=wt, style=st, bbox=bbox)

def pill(x, y, w, h, label, fc, ec="none", sz=7.5, txt_col="white"):
    ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h, boxstyle="round,pad=0.01,rounding_size=0.04",
                 fc=fc, ec=ec, lw=0.8))
    T(x, y, label, sz=sz, ha="center", va="center", col=txt_col, wt="bold")

# ===== Header band =====
box(0, 12.0, 10, 1.0, fc=C_NAVY, ec=C_NAVY, lw=0, r=0.0, style="square,pad=0")
T(0.3, 12.66, "Physics-gated de novo design of VEGF-A miniprotein binders", sz=17, col="white", wt="bold")
T(0.3, 12.18, "A vdW gate exposes ipTM-favoured interpenetration  ·  a 1FLT template rescues a secreted antigen",
  sz=10.5, col="#cfe3f2")

# ===== Pipeline chips =====
chips = ["RFdiffusion\ngeneration", "Developability\nredesign (WALTZ/APR)",
         "MM-GBSA vdW\ngate", "Planned wet-lab\nvalidation"]
cx = [0.35, 2.85, 5.35, 7.85]
for c, x in zip(chips, cx):
    box(x, 11.05, 2.1, 0.78, fc="#e8eef3", ec=C_NAVY, lw=1.0)
    T(x+1.05, 11.44, c, sz=8.6, ha="center", va="center", col=C_NAVY, wt="bold")
for x in cx[:-1]:
    ax.add_patch(FancyArrowPatch((x+2.1, 11.44), (x+2.35, 11.44),
                 arrowstyle="-|>", mutation_scale=10, color=C_NAVY, lw=1.5))

# ===== Finding 1 (left) =====
box(0.15, 6.35, 4.75, 4.45, fc="#ffffff", ec="#333", lw=1.0)
T(0.32, 10.58, "Finding 1  ·  ipTM is fooled by interpenetration", sz=11.5, wt="bold", col=C_V2)
T(0.32, 10.30, "vdW gate catches the artifact in a single step", sz=9, col=C_GREY, st="italic")

# panel (a): ipTM ranking (no gate)  -- x range 0.30..2.40, y range 8.70..10.00
# baseline y=8.70, height up to 9.85 (0..0.90 maps to ipTM 0..0.85)
def bar(ax, x0, x1, base_y, h_max, bars, labels, title, fail_idx=None, y_label=""):
    width_each = (x1-x0)/len(bars)
    for i, (v, c) in enumerate(bars):
        cx = x0 + (i+0.5)*width_each
        bw = width_each*0.45
        h = v/h_max*(h_max-0.15)
        ax.add_patch(Rectangle((cx-bw/2, base_y), bw, h, fc=c, ec="#222", lw=0.6))
        T(cx, base_y+h+0.04, f"{v:.2f}", sz=8, ha="center", wt="bold")
        T(cx, base_y-0.18, labels[i], sz=8, ha="center")
    T((x0+x1)/2, base_y+h_max+0.10, title, sz=9, ha="center", wt="bold", col=C_GREY)
    if fail_idx is not None:
        cx = x0 + (fail_idx+0.5)*width_each
        pill(cx, base_y+0.18, 0.55, 0.18, "FAILS", C_V2)

# bar chart (a) on main ax using data coords
# ipTM ranking: V1 0.78, V2 0.77, V4 0.66 (h_max 0.85)
ip_bars = [(0.78, C_V1), (0.77, C_V2), (0.66, C_V4)]
bar(ax, 0.30, 2.40, 8.65, 0.95, ip_bars, ["V1","V2","V4"], "ipTM (gate OFF)", fail_idx=1)
# V2-looks-best note: short pill (no bbox overlap with chart title), placed
# in the empty space below chart and above the explanatory line below it.
# Sized to fit beside V1/V2/V4 labels without touching the title at y=9.70.
T(1.05, 8.95, "ipTM favours V2", sz=7.0, ha="center", va="center",
  st="italic", col=C_V2,
  bbox=dict(boxstyle="round,pad=0.10", fc="#fde2e4", ec=C_V2, lw=0.7))
ax.annotate("", xy=(1.35, 9.30), xytext=(1.35, 9.05),
            arrowprops=dict(arrowstyle="->", color=C_V2, lw=0.8))

# panel (b): vdW gate -- on main ax x 2.60..4.75
# Map vdW: V1 -43.5 (small), V2 +4803 (clipped), V4 -48.7 (small)
# Use a clipped scale: -50..+50, with V2 shown as full bar in red labelled off-scale
vdw_bars = [(-43.5, C_V1, "-43.5"), (+4803, C_V2, "+4803"), (-48.7, C_V4, "-48.7")]
width_each = (4.75-2.60)/3
for i, (v, c, lab) in enumerate(vdw_bars):
    cx = 2.60 + (i+0.5)*width_each
    bw = width_each*0.45
    h = 0.75  # all same visible height; label carries value
    ax.add_patch(Rectangle((cx-bw/2, 8.65), bw, h, fc=c, ec="#222", lw=0.6))
    T(cx, 8.65+h+0.04, lab, sz=8, ha="center", wt="bold")
    T(cx, 8.65-0.18, ["V1","V2","V4"][i], sz=8, ha="center")
    if i == 1:
        pill(cx, 8.65+0.40, 0.55, 0.18, "FAILS", C_V2)
T(3.675, 8.65+0.75+0.15, "vdW term kcal/mol (gate ON)", sz=9, ha="center", wt="bold", col=C_GREY)

T(0.32, 8.20, "Boltz-1 and OF3 both rank V2's interpenetrating pose as the strongest binder",
  sz=8.4, col="#333")
T(0.32, 7.97, "(ipTM 0.765; vdW +4803 / +4527). The gate rejects it, yielding", sz=8.4, col="#333")
T(0.32, 7.74, "V1  >=  V4  >  V2  across every predictor.", sz=8.6, col=C_NAVY, wt="bold")

# ===== Finding 2 (right) =====
box(5.10, 6.35, 4.75, 4.45, fc="#ffffff", ec="#333", lw=1.0)
T(5.27, 10.58, "Finding 2  ·  template beats MSA", sz=11.5, wt="bold", col=C_V1)
T(5.27, 10.30, "for a secreted, partly-disordered antigen", sz=9, col=C_GREY, st="italic")

# V1 OF3 ipTM bar chart on main ax: x 5.27..9.65
ip2 = [(0.22, C_GREY, "0.22"), (0.76, C_V1, "0.76"), (0.25, C_ACC, "0.25")]
labs2 = ["seq-only", "1FLT tmpl", "ColabFold MSA"]
width_each = (9.65-5.27)/3
for i, (v, c, lab) in enumerate(ip2):
    cx = 5.27 + (i+0.5)*width_each
    bw = width_each*0.5
    h = v/0.85*0.80
    ax.add_patch(Rectangle((cx-bw/2, 8.65), bw, h, fc=c, ec="#222", lw=0.6))
    T(cx, 8.65+h+0.04, lab, sz=8, ha="center", wt="bold")
    # MSA bar shows seed-range inside (room since bar is short and wide)
    if i == 2:
        T(cx, 8.65+h/2, "0.07-0.35", sz=7.5, ha="center", va="center", wt="bold", col="white")
    T(cx, 8.65-0.18, labs2[i], sz=7.6, ha="center")
T(7.46, 8.65+0.80+0.15, "OpenFold3 ipTM of V1", sz=9, ha="center", wt="bold", col=C_GREY)
# ~3.5x arrow from seq bar to tmpl bar
ax.add_patch(FancyArrowPatch((5.27+0.5*width_each+0.10, 8.65+0.22),
             (5.27+1.5*width_each-0.10, 8.65+0.76-0.10),
             arrowstyle="-|>", mutation_scale=12, color=C_V1, lw=1.6))
T(5.27+1.0*width_each+0.05, 8.65+0.55, "~3.5x", sz=9, col=C_V1, wt="bold")
# (replaced by inline label inside MSA bar)

T(5.27, 8.20, "Sequence-only OF3 ipTM for V1 is 0.21-0.22; the 1FLT template lifts it to 0.75-0.77.",
  sz=8.4, col="#333")
T(5.27, 7.97, "ColabFold MSA is seed-unstable (0.07-0.35) and gives no consistent gain over",
  sz=8.4, col="#333")
T(5.27, 7.74, "sequence-only. Template injection is the decisive, model-agnostic fix.", sz=8.4, col="#333")

# ===== Bottom band =====
box(0.15, 0.20, 9.70, 5.95, fc=C_NAVY, ec=C_NAVY, lw=0, r=0.0, style="square,pad=0")
T(0.35, 5.55, "Outcome", sz=13, col="white", wt="bold")
T(0.35, 5.05, "vdW-gated rank:  V1  >=  V4  >  V2", sz=13, col=C_ACC, wt="bold")
T(0.35, 4.55, "V1 first or tied-first in every predictor; V2 rejected by the vdW gate in every mode.",
  sz=9, col="#cfe3f2")
T(0.35, 4.35, "V4 additionally bridges the VEGF dimer", sz=9, col="#cfe3f2")
T(0.35, 4.10, "(chain-pair ipTM 0.66-0.68 vs ~0.07-0.12 for V1/V2).", sz=9, col="#cfe3f2")

# PyMOL-rendered outcome panels (replace the schematic motif)
im1 = imread(os.path.join(_HERE, "GA_V1_contact.png"))
im4 = imread(os.path.join(_HERE, "GA_V4_bridge.png"))
# V1 panel: x 5.55..7.45, y 1.20..3.10  (1.9 x 1.9)
ax.imshow(im1, extent=[5.55, 7.45, 1.20, 3.10], aspect='auto', zorder=3)
ax.add_patch(Rectangle((5.55, 1.20), 1.90, 1.90, fc='none', ec='#5b7fa6', lw=0.8, zorder=4))
T(5.55+0.95, 3.45, "V1 + chain A", sz=8.5, ha="center", va="bottom", wt="bold", col="#cfe3f2")
T(5.55+0.95, 3.26, "weak contact", sz=7.3, ha="center", va="bottom", st="italic", col="#cfe3f2")
T(5.55+0.95, 1.05, "Boltz best seed", sz=7, ha="center", va="top", col="#9bb8d3", st="italic")
# V4 panel: x 7.65..9.55, y 1.20..3.10
ax.imshow(im4, extent=[7.65, 9.55, 1.20, 3.10], aspect='auto', zorder=3)
ax.add_patch(Rectangle((7.65, 1.20), 1.90, 1.90, fc='none', ec='#5b7fa6', lw=0.8, zorder=4))
T(7.65+0.95, 3.45, "V4 bridges VEGF-A dimer", sz=8.5, ha="center", va="bottom", wt="bold", col=C_ACC)
T(7.65+0.95, 3.26, "touches A & B", sz=7.3, ha="center", va="bottom", st="italic", col=C_ACC)
T(7.65+0.95, 1.05, "OF3-template best seed", sz=7, ha="center", va="top", col="#9bb8d3", st="italic")

T(0.35, 0.55, "Reproducible on a single consumer GPU (RTX 5090)  ·  open pipeline, data and code released",
  sz=8.5, col="#9bb8d3", st="italic")

plt.savefig(OUT, dpi=150, facecolor="white")
print("saved", OUT)
