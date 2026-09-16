# -*- coding: utf-8 -*-
"""Supplementary Figure S1. APR-Score featurizer recovery chain and its verification ladder.

Every number printed on the figure is grounded in a disk artefact:
  * 脚本/07_apr_official_featurizer.py        -> recovered scales, 60/60 constraint check
  * 脚本/08_apr_gap_quantification.py         -> official LR coefficients, variance shares, |dP| bounds
  * 脚本/09_ccs_reproduction.py               -> CCS closed form, per-window residuals, terminus conventions
  * 脚本/10_apr_official_integration.py       -> 4 official models repaired, r / |dP| / rank displacement
  * 脚本/11_apr_official_beforeafter.py       -> official-weight nHigh recount on the CRO constructs
  * 脚本/12_sh_entr_ccs_impact.py             -> SH_ENTR sweep (min/mean/max), ordering invariance
  * 脚本/13_ccs_terminus_diag.py              -> N-vs-C terminus asymmetry (3.1e-05 vs 5.8e-02)
  * 脚本/14_ccs_cterm_scan.py                 -> 125 templates + 10 variants, none hit the C-term target
  * 脚本/15_recover_sh_entr.py                -> 24-combination grid, max |r| = 0.5086
  * 脚本/16_sh_entr_granularity.py            -> 1/600 lattice, identifiability verdict
  * tools/aprscore/shc_runtime_evidence.txt   -> byte-identity, strace, shc source comparison

Style: CSBJ/Elsevier-ish, 300 dpi PNG + SVG, English-only text.

PRINT SIZE  (why the layout and the type scale look the way they do)
--------------------------------------------------------------------
This figure is authored 13.6 x 9.6 in and placed 6.5 in wide, so its in-figure
type prints at `fontsize * 6.5 / W_tight`, where W_tight is the tight-cropped
canvas width (13.42 in here, because savefig.bbox is "tight").  The 7 pt floor
therefore demands in-figure type of ~14.9 pt -- the type must be sized for the
PRINTED page, not for the screen.  A 1.16 in x-height at 14.9 pt is the reason
the panels are as sparse as they are.

The first version of this figure used 6.85-8.5 pt type and carried two long
wrapped prose paragraphs.  Lifting it to the floor grew the canvas from 13.8 in
to 28.1 in (2.03x) and made every label collide: those paragraphs were the
artists setting the canvas width, and enlarging them widened the canvas faster
than it widened the type.  Measured, in that order:

    as-is                        canvas 13.817 in -> diverged 2.03x, 98 collisions
    minus the 2 prose paragraphs canvas 10.857 in -> converges, 37 collisions
    + impact rows to one line    -> 10 collisions
    + stage bodies shortened     ->  6 collisions

Both paragraphs now live in the figure caption (Figure S1 in Manuscript_v7.md),
and every remaining label is authored ABOVE the 14.87 pt floor rather than
being lifted to it: lifting flattens the type scale (12.5/13/14/15 all collapse
to 14.87), so the sizes below are chosen so that `figstyle.normalize` finds
nothing to lift.  It is still called before saving as a guard -- it asserts that
nothing prints below 7.2 pt and that no two labels overlap.
"""
import os
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

if len(sys.argv) > 1:
    FIG = sys.argv[1]
else:
    FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(FIG, exist_ok=True)

DPI = 300

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 15.5,
    "figure.dpi": DPI,
    "savefig.dpi": DPI,
    "savefig.bbox": "tight",
})

# palette (matches make_figures.py)
COL_OK    = "#1a6fb5"   # recovered / usable         (blue)
COL_FAIL  = "#c0392b"   # unrecoverable / blocked    (red)
COL_PART  = "#e9a13b"   # partial / bounded          (amber)
COL_NEU   = "#6c757d"   # neutral / context          (grey)
COL_GREEN = "#2a9d8f"   # verified closed            (teal)
FILL_OK   = "#eaf3fb"
FILL_FAIL = "#fdecea"
FILL_PART = "#fdf4e3"
FILL_NEU  = "#f4f5f7"
FILL_GRN  = "#e6f4f2"

# Type scale, in points ON THE AUTHORED CANVAS.  All values sit above the
# 13.98 pt floor so that nothing is lifted and the hierarchy survives.
# 15.0 pt is 1.33 units/char at this figure's scale; card titles are kept to
# <= 17 chars so they fit the 23.2-unit inner width of a card.
FS_LETTER = 19.0   # the A/B/C/D panel letters
FS_TITLE = 17.5    # the four section headings
FS_CARD = 15.0     # stage-card titles (bold + colour adds the emphasis)
FS = 15.0          # everything else


def save(fig, name):
    """Save with the printable-minimum guard applied first."""
    # figstyle.py ships next to this script in BOTH layouts (the authoring copy
    # in 论文/ and the published copy in code/figure_generation/), so the import
    # has to resolve from the script's own directory -- not its parent, which in
    # the authoring tree is the project root and holds no figstyle.py.
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    import figstyle                                     # noqa: E402
    rep = figstyle.normalize(fig, name + ".png", out_dir=FIG, dpi=DPI)
    fig.savefig(os.path.join(FIG, name + ".png"))
    fig.savefig(os.path.join(FIG, name + ".svg"))
    plt.close(fig)
    print("saved", name, "|", figstyle.report_line(rep), flush=True)
    for ov in rep["overlaps_after"]:
        print("   !! text overlap %.2f: %r <> %r" % (ov[3], ov[0], ov[1]),
              flush=True)


def box(ax, x, y, w, h, text, fc=FILL_NEU, ec=COL_NEU, fs=FS, bold=False,
        ls="-", lw=1.1, va="center", ha="center"):
    b = FancyBboxPatch((x, y), w, h,
                       boxstyle="round,pad=0.06,rounding_size=0.10",
                       linewidth=lw, edgecolor=ec, facecolor=fc, linestyle=ls,
                       zorder=2)
    ax.add_patch(b)
    ax.text(x + w / 2 if ha == "center" else x + 0.10, y + h / 2, text,
            ha=ha, va=va, fontsize=fs,
            fontweight="bold" if bold else "normal", zorder=3,
            linespacing=1.35)
    return (x, y, w, h)


def arrow(ax, p1, p2, color=COL_NEU, lw=1.3, ls="-", style="-|>", rad=0.0, ms=9):
    a = FancyArrowPatch(p1, p2, arrowstyle=style, mutation_scale=ms,
                        linewidth=lw, color=color, linestyle=ls, zorder=1,
                        connectionstyle="arc3,rad=%.2f" % rad,
                        shrinkA=1.0, shrinkB=1.0)
    ax.add_patch(a)


def fig_s1():
    fig, ax = plt.subplots(figsize=(13.6, 9.6))
    ax.set_xlim(0, 136)
    ax.set_ylim(0, 92)
    ax.axis("off")

    # ================================================================
    # Panel A (top): the official pipeline, split into the part that
    # runs and the part that cannot run.
    # The "✗ off-host" / "✓ 4/4 load" verdicts are the third line of their own
    # boxes: at 15.5 pt a free-standing label no longer fits between a box top
    # and the section rule.
    # ================================================================
    ax.text(2.0, 84.2, "A", fontsize=FS_LETTER, fontweight="bold", color=COL_NEU, va="top")
    ax.text(6.6, 83.4, "The official APR-Score pipeline, split at its failure point",
            fontsize=FS_TITLE, fontweight="bold", va="top", color="#222222")
    ax.plot([2.0, 134.0], [80.4, 80.4], color="#cccccc", lw=0.8, zorder=1)

    yA = 69.0
    hA = 9.6
    cy = yA + hA / 2                     # every box shares this centre line
    box(ax, 2.0, yA, 12.0, hA, "Sequence\n(6-mer)",
        fc=FILL_NEU, ec=COL_NEU, fs=FS)
    arrow(ax, (14.0, cy), (17.4, cy))

    box(ax, 17.4, yA, 30.0, hA,
        "APR-Score.exe\nSHC-encrypted ELF\nno -r  \u2717 off-host",
        fc=FILL_FAIL, ec=COL_FAIL, fs=FS, lw=1.5)
    arrow(ax, (47.4, cy), (50.8, cy), color=COL_FAIL, lw=1.4)

    box(ax, 50.8, yA, 22.0, hA, "18-D hexapeptide\nfeature vector",
        fc=FILL_PART, ec=COL_PART, fs=FS, lw=1.4)
    arrow(ax, (72.8, cy), (76.2, cy), color=COL_OK, lw=1.4)

    box(ax, 76.2, yA, 27.0, hA,
        "LR / RF / GB / SVM\n.pkcls (Orange)\n\u2713 4/4 load",
        fc=FILL_GRN, ec=COL_GREEN, fs=FS, lw=1.5)
    arrow(ax, (103.2, cy), (106.6, cy), color=COL_OK, lw=1.4)

    box(ax, 106.6, yA, 14.0, hA, "P(AMY)\n= 1",
        fc=FILL_GRN, ec=COL_GREEN, fs=FS, bold=True, lw=1.4)

    # ================================================================
    # Panel B (middle): the recovery chain, 5 stages.
    # Card titles are <= 17 chars and bodies <= 18 chars per line: that is the
    # card's inner width at 15.5 pt.  Every number is retained verbatim; only
    # the framing sentences moved to the caption.
    # ================================================================
    ax.text(2.0, 58.9, "B", fontsize=FS_LETTER, fontweight="bold", color=COL_NEU, va="top")
    ax.text(6.6, 58.1, "The recovery chain, and what each stage closed",
            fontsize=FS_TITLE, fontweight="bold", va="top", color="#222222")
    ax.plot([2.0, 134.0], [55.1, 55.1], color="#cccccc", lw=0.8, zorder=1)

    stages = [
        dict(x=2.0,  title="1  Prop. scales",
             body="4 scales, 20 res.\n60/60 constraints\nmax dev 6.7e-5",
             verdict="CLOSED", vcol=COL_GREEN, fc=FILL_GRN, ec=COL_GREEN),
        dict(x=28.4, title="2  Aggregates",
             body="closed form from\nwindow-difference\nconstraints",
             verdict="CLOSED", vcol=COL_GREEN, fc=FILL_GRN, ec=COL_GREEN),
        dict(x=54.8, title="3  Tripeptide",
             body="CCS = 1 -\nmean(CSS)\n14/15 exact,\nR2 = 0.99898",
             verdict="CLOSED,\nENDS ASIDE",
             vcol=COL_PART, fc=FILL_PART, ec=COL_PART),
        dict(x=81.2, title="4  SH_ENTR",
             body="no closed form:\n324 combos,\nmax |r| 0.51,\nheld 3.250860",
             verdict="NOT\nIDENTIFIABLE",
             vcol=COL_FAIL, fc=FILL_FAIL, ec=COL_FAIL),
        dict(x=107.6, title="5  Official wts.",
             body="4 models repaired\nfor sklearn 1.3+\nmax diff 0",
             verdict="CLOSED", vcol=COL_GREEN, fc=FILL_GRN, ec=COL_GREEN),
    ]
    yB = 27.0
    hB = 27.0
    wB = 25.0          # 25.0 wide + 1.4 gutter = the 26.4 pitch the x's assume
                       # (the longest title, "1  Property scales", is 22.6 units
                       # at 15 pt, so the inner width has to exceed that)
    for i, s in enumerate(stages):
        box(ax, s["x"], yB, wB, hB, "", fc=s["fc"], ec=s["ec"], lw=1.4)
        ax.text(s["x"] + 1.2, yB + hB - 1.4, s["title"], fontsize=FS_CARD,
                fontweight="bold", va="top", ha="left", color=s["ec"])
        ax.text(s["x"] + 1.2, yB + hB - 4.8, s["body"], fontsize=FS,
                va="top", ha="left", color="#333333", linespacing=1.35)
        # verdict chip
        ax.add_patch(Rectangle((s["x"] + 1.2, yB + 1.0), wB - 2.4, 6.0,
                               facecolor=s["vcol"], edgecolor="none", alpha=0.16, zorder=2))
        ax.text(s["x"] + wB / 2, yB + 4.0, s["verdict"], fontsize=FS,
                fontweight="bold", ha="center", va="center", color=s["vcol"],
                linespacing=1.25, zorder=3)
        if i < len(stages) - 1:
            arrow(ax, (s["x"] + wB, yB + hB / 2), (stages[i + 1]["x"], yB + hB / 2),
                  color=COL_NEU, lw=1.2, ms=8)

    # ================================================================
    # Panel C (bottom-left): verification ladder.
    # Row pitch 3.0 units == one line of 15.5 pt type plus leading.
    # ================================================================
    ax.text(2.0, 25.5, "C", fontsize=FS_LETTER, fontweight="bold", color=COL_NEU, va="top")
    ax.text(6.6, 24.7, "Verification ladder",
            fontsize=FS_TITLE, fontweight="bold", va="top", color="#222222")
    ax.plot([2.0, 60.0], [22.5, 22.5], color="#cccccc", lw=0.8, zorder=1)

    ladder = [
        ("1", "60/60 constraints", "6.7 \u00d7 10\u207b\u2075", COL_GREEN),
        ("2", "LR hand vs sklearn", "0.00 \u00d7 10\u2070", COL_GREEN),
        ("3", "CCS per window 14/15", "5 \u00d7 10\u207b\u2075", COL_GREEN),
        ("4", "CCS N-terminus", "3.1 \u00d7 10\u207b\u2075", COL_GREEN),
        ("5", "CCS C-terminus", "5.8 \u00d7 10\u207b\u00b2  \u2717", COL_FAIL),
        ("6", "SH_ENTR 324 combos", "|r| = 0.51  \u2717", COL_FAIL),
        ("7", "Order of 6 constructs", "invariant", COL_OK),
    ]
    yC = 20.3
    for i, (n, what, res, col) in enumerate(ladder):
        yy = yC - i * 3.0
        # 2.4 / 6.0, not 2.6 / 4.9: at 15.5 pt the "n." glyph box is ~2.1 units
        # wide and overlapped the label by 13 % of its own area at 4.9.
        ax.text(2.4, yy, n + ".", fontsize=FS, fontweight="bold", color=col, va="center")
        ax.text(6.0, yy, what, fontsize=FS, color="#333333", va="center")
        ax.text(60.0, yy, res, fontsize=FS, color=col, va="center",
                ha="right", fontweight="bold")

    # ================================================================
    # Panel D (bottom-right): impact, bounded.
    # One line per row -- label (with its status folded in) plus a short detail.
    # The former two-line label/tag stack plus a three-line detail cannot fit
    # four rows at 15.5 pt; the dropped qualifiers are in the caption.
    # ================================================================
    ax.text(64.0, 25.5, "D", fontsize=FS_LETTER, fontweight="bold", color=COL_NEU, va="top")
    ax.text(68.6, 24.7, "What the residual gap does to the conclusions",
            fontsize=FS_TITLE, fontweight="bold", va="top", color="#222222")
    ax.plot([64.0, 134.0], [22.5, 22.5], color="#cccccc", lw=0.8, zorder=1)

    # Label = what it is (left, bold); detail = status + numbers (right).
    # Each label is <= 22 chars and each detail <= 33 chars, which at 15.5 pt
    # leaves >= 4 units between the two columns at x = 66 and x = 134.
    impact = [
        ("Score ordering", "invariant; r = 0.986, rho 0.97", COL_GREEN),
        ("Window counts", "bounded: V1 0-1, V2 4-20, V4 0-9", COL_PART),
        ("Reported nHigh", "official wts. on 17-D", COL_PART),
        ("CCS at the C-terminus", "~0.004 in P; below 0.5", COL_NEU),
    ]
    yD = 19.6
    for i, (what, detail, col) in enumerate(impact):
        yy = yD - i * 3.4
        ax.add_patch(Rectangle((64.4, yy - 1.4), 0.6, 2.8,
                               facecolor=col, edgecolor="none", zorder=2))
        ax.text(66.0, yy, what, fontsize=FS, fontweight="bold",
                color="#222222", va="center")
        ax.text(134.0, yy, detail, fontsize=FS, color="#444444",
                va="center", ha="right")

    fig.tight_layout()
    save(fig, "FigS1_apr_recovery")


if __name__ == "__main__":
    fig_s1()
