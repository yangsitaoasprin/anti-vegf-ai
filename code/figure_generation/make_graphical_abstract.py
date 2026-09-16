# =============================================================================
# Author : Yang Sitao (杨四涛)
# Affiliation : School of Pharmaceutical Sciences, Dali University;
#               The Third Affiliated Hospital of Dali University
# Contact : yangsitaoasprin@swpu.edu.cn
# License : MIT (see ../../LICENSE)
# Purpose: Assembles the Graphical Abstract (matplotlib) using the two PyMOL
#          panels produced by render_ga_outcome.py.
# =============================================================================
"""Graphical Abstract for CSBJ submission.

WHY IT LOOKS LIKE THIS
----------------------
Elsevier's graphical-abstract guidance (elsevier.com/researcher/author/
tools-and-resources/graphical-abstract) requires, in its own words:

  * "a minimum of 1328 x 531 pixels (w x h) using a minimum resolution of
     300 dpi. If you are submitting a larger image, please use the same ratio
     (500 wide x 200 high)"
  * "your image will be scaled proportionally to fit in the available window
    on ScienceDirect: a 500 by 200-pixel rectangle"
  * "Font: please use Times, Arial, Courier or Symbol with a large enough font
     size as the image will be reduced in size for the table of contents to fit
     a window 200 pixels high"

So the binding constraints are (a) a 2.5:1 LANDSCAPE frame and (b) legibility
after reduction to a 200-px-high window -- not the 7 pt print floor that governs
figures in the article body.  The previous version was a 10 x 13 in PORTRAIT
poster with ~30 lines of prose; squeezed into a 500 x 200 window it would have
been reduced to ~154 x 200 px and been unreadable.

The frame here is 5 x 2 in = 1500 x 600 px at 300 dpi: the required 2.5:1 ratio,
above the 1328 x 531 minimum.  Every label is >= 9 pt, i.e. >= 12.5 px once the
width is reduced to Elsevier's 500-px window, and ~11.2 pt at the 6.2 in width
build_docx.py places it at (FIG_WIDTH["GA"]) -- above the 7 pt figure floor,
which the GA is not required to meet but which makes the same file safe to reuse
in the article.  figstyle reports a slightly conservative 10.7 pt because it
measures the tight-bbox width (5.20 in) while this file is saved at exactly 5.0
in.  Layout is one 0..10 x 0..4 grid with 1 grid unit = 0.5 in, so text sizes
and box positions are directly comparable.
"""
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.image import imread

_HERE = os.path.dirname(os.path.abspath(__file__))


def _find_paper_dir():
    """Manuscript folder (<project>/论文), located by walking up from here.

    Two layouts must work with this same file:
      * authoring tree -- .../repro/code/figure_generation/, with the figures
        consumed by build_docx.py living in <project>/论文/figures;
      * standalone clone of the published repro repository -- code/ sits at the
        repo root and there is no 论文/ folder at all, so this returns None and
        the caller falls back to a figures/ dir inside the clone.  Writing
        outside the clone (the old dirname x3 arithmetic) is never correct
        there: it lands in the clone's parent directory.
    """
    here = _HERE
    for _ in range(6):
        if os.path.isdir(os.path.join(here, "论文")):
            return os.path.join(here, "论文")
        parent = os.path.dirname(here)
        if parent == here:
            break
        here = parent
    return None


PAPER_DIR = _find_paper_dir()
if PAPER_DIR is None:
    FIGURES_DIR = os.path.join(os.path.dirname(os.path.dirname(_HERE)),
                               "figures")
else:
    FIGURES_DIR = os.path.join(PAPER_DIR, "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)
OUT = os.path.join(FIGURES_DIR, "Graphical_Abstract.png")

C_V1 = "#2a6f97"; C_V2 = "#c1121f"; C_V4 = "#2a9d8f"
C_ACC = "#e9c46a"; C_GREY = "#6c757d"; C_NAVY = "#0d3a5c"

# Arial / Helvetica / Liberation Sans first: Elsevier's GA spec lists only
# Times, Arial, Courier and Symbol, and a bare sans-serif resolves to whatever
# the build machine has.  DejaVu is the guaranteed matplotlib fallback.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "Liberation Sans",
                                   "DejaVu Sans"]

W_IN, H_IN = 5.0, 2.0                 # 2500 x 1000 px @ 500 dpi -> 2.5 : 1
# Elsevier GA spec: 1328x531 px @ 300 dpi minimum.  At placed 6.2 in we need
# at least 6.2*300 = 1860 px on the short edge; 500 dpi gives 2500 px and
# leaves headroom against the 1328x531 floor.
DPI = 500
fig = plt.figure(figsize=(W_IN, H_IN), dpi=DPI)
ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
ax.set_xlim(0, 10); ax.set_ylim(0, 4)

# ---------------------------------------------------------------- primitives
def T(x, y, s, sz=9.5, col="#1a1a1a", ha="left", va="bottom", wt="normal",
      st="normal", zorder=5):
    return ax.text(x, y, s, fontsize=sz, color=col, ha=ha, va=va, weight=wt,
                   style=st, zorder=zorder)

def box(x, y, w, h, fc, ec="none", lw=1.0, r=0.05):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                 boxstyle="round,pad=0,rounding_size=%.3f" % r,
                 fc=fc, ec=ec, lw=lw, zorder=1))

def pill(x, y, w, h, label, fc, sz=9.0):
    ax.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h,
                 boxstyle="round,pad=0,rounding_size=0.05",
                 fc=fc, ec="none", zorder=6))
    # zorder 7: T() defaults to 5, i.e. BELOW the patch added on the line above,
    # so the pill used to paint over its own white label and the label vanished.
    # `w` must also exceed the label's own width (see the call site).
    T(x, y, label, sz=sz, ha="center", va="center", col="white", wt="bold",
      zorder=7)

def card(x0, x1, heading, hcol):
    """White evidence card with a one-line heading."""
    box(x0, 1.30, x1 - x0, 1.92, fc="#ffffff", ec="#b8c4cc", lw=0.9, r=0.07)
    T(x0 + 0.14, 2.90, heading, sz=9.0, wt="bold", col=hcol)

def minibars(x0, x1, vals, colors, ticklabels, vlabels, vmax, base=1.74,
             hmax=0.58):
    """Bars up from `base`; value above each bar, tick label below."""
    n = len(vals); slot = (x1 - x0) / n
    bw = slot * 0.44
    for i, v in enumerate(vals):
        cx = x0 + (i + 0.5) * slot
        h = (v / vmax) * hmax
        ax.add_patch(Rectangle((cx - bw/2, base), bw, h, fc=colors[i],
                               ec="#222", lw=0.6, zorder=3))
        T(cx, base + h + 0.05, vlabels[i], sz=9.0, ha="center", wt="bold",
          col="#222", zorder=4)
        T(cx, base - 0.06, ticklabels[i], sz=9.0, ha="center", va="top",
          col="#222", zorder=4)
    return slot

# ============================================================ header band
box(0, 3.30, 10, 0.70, fc=C_NAVY, r=0)
T(0.20, 3.62, "Physics-gated AI miniprotein design on VEGF-A",
  sz=13.0, col="white", wt="bold")
T(0.20, 3.36, "vdW gate beats ipTM  \u00b7  1FLT template beats MSA",
  sz=9.5, col="#cfe3f2")

# ===================================================== card 1: ipTM ranking
card(0.16, 3.28, "1 \u00b7 ipTM picks V2", C_V2)
# Boltz-1 full-length 5-seed mean ipTM: V1 0.464, V2 0.521, V4 0.357
minibars(0.34, 3.10, [0.464, 0.521, 0.357],
         [C_V1, C_V2, C_V4], ["V1", "V2", "V4"],
         ["0.46", "0.52", "0.36"], vmax=0.95)

# ============================================== card 2: the vdW gate (physics)
card(3.42, 6.54, "2 \u00b7 vdW gate", C_GREY)
# Signed scale clipped at +/-1 drawn unit: V2's +4803 is off-scale, so it is
# drawn clipped with its true value printed above (as in Fig. 6).
_gx0, _gx1 = 3.60, 6.36
_slot = (_gx1 - _gx0) / 3
_zero = 1.92                       # y of the vdW = 0 line
# Drawn height per 50 kcal/mol of |vdW|.  0.24 keeps the tallest NEGATIVE bar
# (-48.7) clear of the tick-label band at y = 1.45..1.68 -- the band cards 1 and
# 3 use.  This is the only card whose bars hang below their own zero line, so
# the unit size has to be chosen against that band.  V2's bar is clipped at a
# fixed height and is unaffected.
_ax_unit = 0.24
for _i, (_v, _c, _lab, _tl) in enumerate([
        (-43.5, C_V1, "-43.5", "V1"),
        (+4803, C_V2, "+4803", "V2"),
        (-48.7, C_V4, "-48.7", "V4")]):
    _cx = _gx0 + (_i + 0.5) * _slot
    _bw = _slot * 0.44
    if _v > 0:
        _h = _ax_unit * 1.05                      # clipped at the top
        ax.add_patch(Rectangle((_cx - _bw/2, _zero), _bw, _h, fc=C_V2,
                               ec="#222", lw=0.6, hatch="///", zorder=3))
        # 1.50 x 0.30 units: "REJECTED" measures 1.34 x 0.23 units at 9 pt, so
        # the previous 0.78 x 0.24 pill was narrower AND shorter than its label.
        pill(_cx, _zero + _h + 0.22, 1.50, 0.30, "REJECTED", C_V2, sz=9.0)
    else:
        _h = abs(_v) / 50.0 * _ax_unit
        ax.add_patch(Rectangle((_cx - _bw/2, _zero - _h), _bw, _h, fc=_c,
                               ec="#222", lw=0.6, zorder=3))
    # The clipped bar's true value goes in the card's top band, above the
    # REJECTED pill; the negative bars' values sit just above the zero line.
    T(_cx, 2.56 if _v > 0 else _zero + 0.06, _lab, sz=9.0, ha="center",
      va="bottom", wt="bold", col="#222", zorder=4)
    # Same tick band as cards 1 and 3.  These used to sit 0.38 units lower,
    # i.e. below the card's own bottom edge and on top of the outcome band's
    # headline ("vdW-gated ranking: ..."), which is what the gate flagged.
    T(_cx, 1.68, _tl, sz=9.0, ha="center", va="top", col="#222")
ax.plot([_gx0, _gx1], [_zero, _zero], color="#444", lw=0.9, zorder=2)
# The unit lives in the card's top band, right-aligned.  As a free-standing
# caption at the foot of the card it ran into the -48.7 value label.
T(6.40, 2.90, "(kcal/mol)", sz=9.0, ha="right", col="#666")

# ==================================== card 3: 1FLT template vs ColabFold MSA
card(6.68, 9.80, "3 \u00b7 1FLT vs MSA", C_V1)
# OpenFold3 ipTM of V1: sequence-only 0.22, 1FLT template 0.76, ColabFold MSA 0.25
_s3 = minibars(6.90, 9.58, [0.22, 0.76, 0.25],
               [C_GREY, C_V1, C_ACC], ["seq-only", "1FLT", "MSA"],
               ["0.22", "0.76", "0.25"], vmax=0.95)
# From just above the grey seq-only bar to the top-left corner of the 1FLT bar.
# It used to end at (x + 0.30, 1.74 + 0.76 - 0.06) = (8.54, 2.44), i.e. inside
# the "0.76" value label, which the shaft was drawn straight through.  The text
# gate cannot see this: it compares text with text, and the arrow is a patch.
ax.add_patch(FancyArrowPatch((6.90 + 0.5*_s3 - 0.10, 1.74 + 0.16),
                             (6.90 + 1.5*_s3 - 0.18, 1.74 + 0.46),
                             arrowstyle="-|>", mutation_scale=9,
                             color=C_V1, lw=1.4, zorder=5))
# Raised clear of the 1FLT bar's own value label ("0.76", which tops out at
# y = 2.48); at 0.60 it was printed straight through it.
T(6.90 + 1.5*_s3, 1.74 + 0.80, "\u00d73.5", sz=9.5, col=C_V1, wt="bold",
  ha="left", zorder=6)

# ============================================================ outcome band
# The two PyMOL outcome thumbnails (V1 contact / V4 dimer bridge) used to sit
# in the right half of this band at +30 % size; the author decided the figure
# is cleaner without them and asked for them to be removed.  The band stays
# widened (1.10 units, y = 0.30..1.40) to host the three summary lines.
box(0, 0.30, 10, 1.10, fc=C_NAVY, r=0)
# The yellow "vdW-gated ranking..." line used to print at y = 1.18 (baseline),
# whose glyph top then crossed into the card band (y >= 1.30) and overpainted
# the V1 / V2 / V4 tick labels below the three bar charts.  Pulling the
# baseline to y = 1.14 keeps the glyph top at 1.30 -- flush with the card
# ceiling and clear of the tick labels.
T(0.20, 1.14, "vdW-gated ranking:  V1  \u2265  V4  >  V2",
  sz=11.5, col=C_ACC, wt="bold")
T(0.20, 0.86, "V1 leads in every predictor; V2 fails the gate in every mode.",
  sz=9.5, col="#cfe3f2")
T(0.20, 0.58, "V4 bridges the VEGF-A dimer (ipTM 0.66\u20130.68); V1 does not.",
  sz=9.5, col="#cfe3f2")

# ------------------------------------------------------------------- gate
# figstyle.py ships next to this script in both layouts, so the import does not
# depend on the manuscript folder existing (it does not, in a standalone clone).
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
import figstyle                                        # noqa: E402
_rep = figstyle.normalize(fig, "Graphical_Abstract.png", out_dir=FIGURES_DIR,
                          dpi=DPI, min_pt=9.0)

fig.savefig(OUT, dpi=DPI, facecolor="white")
fig.savefig(os.path.join(FIGURES_DIR, "Graphical_Abstract.pdf"),
            facecolor="white")                          # vector, preferred
fig.savefig(os.path.join(FIGURES_DIR, "Graphical_Abstract.tiff"),
            dpi=DPI, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
print("saved", OUT)
print("   ", figstyle.report_line(_rep))
for _ov in _rep["overlaps_after"]:
    print("   !! overlap %.2f: %r <> %r" % (_ov[3], _ov[0], _ov[1]))
