# -*- coding: utf-8 -*-
# =============================================================================
# Author : Yang Sitao (杨四涛)
# Affiliation : School of Pharmaceutical Sciences, Dali University;
#               The Third Affiliated Hospital of Dali University
# Contact : yangsitaoasprin@swpu.edu.cn
# License : MIT (see ../../LICENSE)
# Purpose: Regenerates Figs. 1-8 of the manuscript from hardcoded real values.
# =============================================================================
"""Publication-quality figures, all data from real results/ outputs.
Style: CSBJ/Elsevier-ish, 300 dpi PNG + SVG, English-only text.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

if len(sys.argv) > 1:
    FIG = sys.argv[1]
else:
    FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

COL_FINAL = "#1a6fb5"
COL_ORIG  = "#c0c0c0"
COL_RED   = "#c0392b"
COL_ACC   = "#2a9d8f"
COL_WARN  = "#e9c46a"

def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".png"))
    fig.savefig(os.path.join(FIG, name + ".svg"))
    plt.close(fig)
    print("saved", name, flush=True)

# =====================================================================
# Figure 1. Pipeline schematic
# =====================================================================
def fig1_pipeline():
    fig, ax = plt.subplots(figsize=(11.5, 6.0))
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

    def box(x, y, w, h, text, fc="#eaf3fb", ec=COL_FINAL, fs=9.5, bold=False):
        b = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08,rounding_size=0.12",
                           linewidth=1.2, edgecolor=ec, facecolor=fc)
        ax.add_patch(b)
        ax.text(x+w/2, y+h/2, text, ha="center", va="center",
                fontsize=fs, fontweight="bold" if bold else "normal", wrap=True)

    def arrow(x1, y1, x2, y2, color="#555"):
        a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=16,
                            linewidth=1.5, color=color)
        ax.add_patch(a)

    def gate_badge(x, y, w, text, color):
        """Coloured 'guard-gate' tag bar sitting on top of a box."""
        b = FancyBboxPatch((x, y), w, 0.30, boxstyle="round,pad=0.03,rounding_size=0.08",
                           linewidth=0, edgecolor=color, facecolor=color)
        ax.add_patch(b)
        ax.text(x+w/2, y+0.15, text, ha="center", va="center", fontsize=7.8,
                fontweight="bold", color="white")

    # Stage labels (high enough above the topmost boxes)
    ax.text(1.35, 9.75, "Stage 1: Generation", fontsize=10.5, fontweight="bold", color=COL_FINAL, ha="center")
    ax.text(4.25, 9.75, "Stage 2: Developability loop", fontsize=10.5, fontweight="bold", color=COL_ACC, ha="center")
    ax.text(7.1, 9.75, "Stage 3: Cross-model validation", fontsize=10.5, fontweight="bold", color="#7b2cbf", ha="center")
    ax.text(9.3, 9.75, "Stage 4", fontsize=10.5, fontweight="bold", color="#2e7d32", ha="center")

    box(0.3, 8.2, 2.1, 1.3, "RFdiffusion\nhotspot-guided binder growth\non VEGF-A RBD (60/80/100 aa)", fc="#eaf3fb")
    box(0.3, 6.2, 2.1, 1.3, "ProteinMPNN\nstructure -> sequence\n(3 temperatures)", fc="#eaf3fb")
    arrow(1.35, 8.2, 1.35, 7.5)
    box(0.3, 4.2, 2.1, 1.3, "ESMFold coarse filter\npLDDT >= 70 & RMSD <= 2 A\nfold-plausibility only", fc="#eaf3fb")
    arrow(1.35, 6.2, 1.35, 5.5)
    box(0.3, 2.2, 2.1, 1.3, "Boltz-1 re-ranking\n105 complexes\nipTM / pAE + 3-seed stability", fc="#eaf3fb")
    arrow(1.35, 4.2, 1.35, 3.5)
    box(0.3, 0.2, 2.1, 1.3, "Composite score\n0.40*I + 0.25*dG + 0.20*pLDDT\n+ 0.15*developability -> Top-6", fc="#eaf3fb")
    arrow(1.35, 2.2, 1.35, 1.5)

    box(3.0, 4.9, 2.5, 1.4, "Developability screen\nWALTZ + APR-Score\namyloid aggregation check", fc="#fff8e1", ec=COL_WARN, bold=True)
    arrow(2.4, 3.0, 3.0, 4.9, color="#bbb")

    box(3.0, 2.5, 2.5, 1.4, "Iterative redesign\nProteinMPNN with fixed\ninterface residues", fc="#e8f6ef", ec=COL_ACC)
    arrow(3.0, 4.9, 3.0, 3.9)
    arrow(3.0, 2.5, 3.0, 1.9, color=COL_ACC)
    ax.text(3.0, 1.8, "re-screen with WALTZ/APR", fontsize=8.5, color=COL_ACC, ha="center", style="italic")

    box(3.0, 0.2, 2.5, 1.3, "Full-length VEGF165 recheck\n(HBD 99-165 context)\n3-seed Boltz + MM-GBSA", fc="#eaf3fb")
    arrow(3.0, 1.9, 3.0, 1.5, color=COL_ACC)
    arrow(5.5, 0.9, 6.0, 0.9, color="#888")  # full-length recheck -> vdW gate

    box(6.0, 4.9, 2.2, 1.4, "OpenFold3 cross-validation\n1FLT template > ColabFold MSA\n(dimer + full-length)", fc="#f3e8fb", ec="#7b2cbf")
    gate_badge(6.0, 6.30, 2.2, "GATE 2 \u2014 template beats MSA", "#7b2cbf")
    arrow(5.5, 5.5, 6.0, 5.5)
    box(6.0, 2.5, 2.2, 1.4, "MM-GBSA vdW gate\n(physical arbiter)\nvdW<0 & dG<0 = clean binding", fc="#fdecea", ec=COL_RED, bold=True)
    gate_badge(6.0, 3.90, 2.2, "GATE 1 \u2014 physics beats ipTM", COL_RED)
    arrow(6.0, 4.9, 6.0, 4.22)
    box(6.0, 0.2, 2.2, 1.3, "vdW-gated ranking\nV1 >= V4 > V2\n(Boltz + MM-GBSA; V2 excluded)", fc="#fef9e7", ec=COL_WARN, bold=True)
    arrow(6.0, 2.5, 6.0, 1.5)

    box(8.6, 0.2, 1.3, 5.9, "Wet lab\n(planned)\n\nfull-length\nVEGF165\nSPR/BLI KD\n\ndual-receptor\nblocking\n\nHUVEC assay", fc="#eef7f0", ec="#2e7d32")
    arrow(8.2, 0.9, 8.6, 0.9)
    save(fig, "Fig1_pipeline")

# =====================================================================
# Figure 2. Redesign improves binding (98-aa dG: original vs final)
# =====================================================================
def fig2_redesign_dG():
    cands = ["V1  #1_T0.3_s132", "V2  #2c_T0.1_s18", "V4  T0.1_s35"]
    orig = np.array([-89.7, -57.3, -37.8])
    final = np.array([-129.9, -53.3, -51.4])
    err_f = np.array([2.2, 17.2, 6.1])
    x = np.arange(len(cands)); w = 0.36

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.bar(x - w/2, orig, w, label="original design (98-aa)", color=COL_ORIG, edgecolor="#999", linewidth=0.7)
    ax.bar(x + w/2, final, w, yerr=err_f, capsize=3, label="final redesigned (98-aa)",
           color=COL_FINAL, edgecolor="#0d3a5c", linewidth=0.7)
    for xi, v in zip(x + w/2, final):
        ax.annotate(f"{v:.0f}", (xi, v + 3), ha="center", fontsize=8.5, fontweight="bold")
    for xi, v in zip(x - w/2, orig):
        ax.annotate(f"{v:.0f}", (xi, v + 3), ha="center", fontsize=8.5, color="#555")
    ax.axhline(0, color="k", linewidth=0.8)
    ax.set_xticks(x); ax.set_xticklabels(cands)
    ax.set_ylabel("MM-GBSA dG (kcal/mol)")
    ax.set_title("Redesign: predicted binding improved (V1, V4) or maintained (V2)\n(more negative = stronger)")
    ax.legend(frameon=False, fontsize=9)
    ax.set_ylim(-150, 5)
    ax.text(0.02, 0.95, "all vdW < 0 (no clash artifacts); V2 = conservative interface fine-tuning\nfor dual-receptor (KDR+Flt-1) blocking, not for ΔG gain",
            transform=ax.transAxes, fontsize=8, color=COL_ACC, style="italic")
    save(fig, "Fig2_redesign_dG")

# =====================================================================
# Figure 3. Full-length recheck (98-aa truncated vs 165-aa full-length)
# =====================================================================
def fig3_full_length():
    cands = ["V1", "V2", "V4"]
    dg98 = np.array([-129.9, -53.3, -51.4])
    dg165 = np.array([-100.8, -55.0, -95.0])   # std-sequence recheck (2026-08-18)
    err98 = np.array([2.2, 17.2, 6.1])
    err165 = np.array([20.6, 32.9, 24.6])      # sample SD of clean seeds (n-1)
    clean165 = ["3/3", "2/3", "2/3"]
    x = np.arange(len(cands)); w = 0.36

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.bar(x - w/2, dg98, w, yerr=err98, capsize=3, label="VEGF-A RBD (1-98)",
           color="#9db8d2", edgecolor="#567", linewidth=0.7)
    ax.bar(x + w/2, dg165, w, yerr=err165, capsize=3, label="full-length VEGF-A165 (with HBD 99-165)",
           color=COL_FINAL, edgecolor="#0d3a5c", linewidth=0.7)
    for xi, v, c in zip(x + w/2, dg165, clean165):
        ax.annotate(f"{v:.0f} ({c})", (xi, v + 3), ha="center", fontsize=8.5, fontweight="bold")
    for xi, v in zip(x - w/2, dg98):
        ax.annotate(f"{v:.0f}", (xi, v + 3), ha="center", fontsize=8.5, color="#555")
    ax.set_xticks(x); ax.set_xticklabels(cands)
    ax.set_ylabel("MM-GBSA dG (kcal/mol)")
    ax.set_title("Full-length recheck: binding remains strong and vdW-clean")
    ax.legend(frameon=False, fontsize=9)
    ax.set_ylim(-150, 5)
    ax.text(0.02, 0.95, "relative order V1 >= V4 > V2; clash seeds (vdW>0) excluded",
            transform=ax.transAxes, fontsize=8.5, color=COL_ACC, style="italic")
    ax.text(0.02, 0.88, "absolute dG not directly comparable between constructs (Section 3.5)",
            transform=ax.transAxes, fontsize=7.5, color="#666", style="italic")
    save(fig, "Fig3_full_length")

# =====================================================================
# Figure 4. Cross-validation dG + per-seed vdW gating
# =====================================================================
def fig4_crossval():
    cands = ["V1", "V4", "V2"]
    boltz = np.array([-100.8, -95.0, -55.0])   # std-sequence full-length recheck
    of3ss = np.array([-121.9, -53.9, -44.2])    # sequence-only (retained)
    of3tmpl = np.array([-115.0, -51.1, 4613.1]) # std-sequence template recheck
    x = np.arange(len(cands)); w = 0.26

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.2), gridspec_kw={"width_ratios": [1, 1.3]})

    # --- Panel (a): dG. V2 OF3-template value (+4613) is off-scale; drawn as hatched
    # placeholder + red ↑ marker so the on-scale bars and the off-scale annotation do
    # not overlap visually.
    ax1.bar(x - w, boltz, w, label="Boltz-1 (full-length, clean-seed mean)", color="#5fa8d3", edgecolor="#246")
    ax1.bar(x, of3ss, w, label="OF3 sequence-only (initial batch)", color="#84a59d", edgecolor="#365")
    ax1.bar(x[:-1] + w, of3tmpl[:-1], w, label="OF3 template (1FLT)", color=COL_FINAL, edgecolor="#0d3a5c")
    ax1.bar(x[2] + w, 240, w, color="none", edgecolor=COL_RED, hatch="///", linewidth=1.0)
    ax1.scatter([x[2] + w], [240], marker="^", s=80, color=COL_RED, zorder=5)
    ax1.axhline(0, color="k", linewidth=0.8)
    ax1.set_xticks(x); ax1.set_xticklabels(cands)
    ax1.set_ylabel("MM-GBSA dG (kcal/mol)")
    ax1.set_title("(a) Binding free energy")
    ax1.legend(frameon=False, fontsize=7.0, loc="lower left")
    ax1.set_ylim(-160, 320)
    ax1.annotate("dG = +4613.1 (off-scale)\nvdW = +4527.0\nsevere clash",
                 (x[2]+w, 240), xytext=(x[2]+w, 300), ha="center", va="top",
                 fontsize=7.5, color=COL_RED, fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color=COL_RED, lw=1.0))

    # --- Panel (b): per-seed vdW scatter.
    # The two interpenetrating Boltz/OF3-template V2 seeds (vdW +4803 / +4527) exceed the
    # panel range; they are drawn as red markers pinned to the panel top with explicit values.
    seeds = {
        "Boltz-1 (3 seeds each)": [(-5.0, 0), (-72.9, 1), (-52.4, 2),
                                   (-33.8, 3), (76.0, 4), (-63.6, 5),
                                   (-0.8, 6), (-50.5, 7)],
        "OF3 sequence-only": [(-72.1, 0), (-59.1, 3), (10.3, 6)],
        "OF3 template (1FLT)": [(-74.8, 0), (-49.4, 3)],
    }
    colormap = {"Boltz-1 (3 seeds each)": "#5fa8d3", "OF3 sequence-only": "#84a59d", "OF3 template (1FLT)": COL_FINAL}
    for method, pts in seeds.items():
        xs, ys, cl = [], [], []
        for v, i in pts:
            xs.append(i); ys.append(v)
            cl.append(COL_RED if v > 0 else colormap[method])
        ax2.scatter(xs, ys, s=60, color=cl, edgecolor="white", linewidth=0.5, zorder=3,
                    label=method)
    # off-scale interpenetrating seeds pinned at the top edge of the panel
    for v, i in [(4803.0, 8), (4527.0, 6)]:
        ax2.scatter([i], [195], s=80, marker="v", color=COL_RED, edgecolor="white",
                    linewidth=0.8, zorder=5)
        ax2.annotate(f"+{v:.0f}", (i, 215), ha="center", va="bottom", fontsize=8,
                     color=COL_RED, fontweight="bold")
    ax2.axhline(0, color="k", linewidth=1.2)
    ax2.axhspan(-160, 0, color=COL_ACC, alpha=0.08)
    for xi in (2.5, 5.5):
        ax2.axvline(xi, color="#bbb", linestyle=":", linewidth=0.8)
    ax2.set_xticks([1, 4, 7]); ax2.set_xticklabels(cands)
    ax2.set_xlim(-0.6, 8.6)
    ax2.set_ylabel("MM-GBSA vdW term (kcal/mol)")
    ax2.set_title("(b) Per-seed vdW gating (vdW > 0 = interpenetrating pose)")
    ax2.set_ylim(-160, 260)
    ax2.legend(frameon=False, fontsize=7.5, loc="upper left")
    ax2.text(0.02, 0.55, "green = clean-binding zone; red triangles = off-scale clash seeds",
             transform=ax2.transAxes, fontsize=8, color=COL_ACC)
    fig.suptitle("Cross-method validation and physical adjudication (MM-GBSA vdW gating)",
                 fontsize=11.5)
    save(fig, "Fig4_crossval")

# =====================================================================
# Figure 5. Developability improvements
# =====================================================================
def fig5_developability():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.0))

    regions = ["TALAIA\n(15-20)", "AIALFAAA\n(64-71)"]
    p_orig = np.array([0.584, 0.548])
    p_new = np.array([0.482, 0.488])
    x = np.arange(2); w = 0.36
    ax1.bar(x - w/2, p_orig, w, label="#4 original (vegf_len80_3)", color=COL_RED, alpha=0.75)
    ax1.bar(x + w/2, p_new, w, label="V4 (16 mutations)", color=COL_ACC)
    for xi, (a, b) in enumerate(zip(p_orig, p_new)):
        ax1.annotate(f"-{100*(a-b)/a:.0f}%", (xi, max(a, b) + 0.02), ha="center",
                     fontsize=9, color="#333", fontweight="bold")
    ax1.set_xticks(x); ax1.set_xticklabels(regions)
    ax1.set_ylabel("APR-Score P(AMY=1)")
    ax1.set_ylim(0.4, 0.68)
    ax1.set_title("(a) APR-Score P(AMY=1) at the WALTZ-positive regions of #4")
    ax1.legend(frameon=False, fontsize=8, loc="lower left")
    ax1.text(0.55, 0.95, "WALTZ high-risk regions removed;\nAPR agrees directionally",
             transform=ax1.transAxes, fontsize=7.5, color="#555", va="top")

    groups = ["#1 original", "V1\n(de-Ala, 16 mutations)"]
    nhigh = np.array([31, 17])
    bars = ax2.bar(groups, nhigh, color=[COL_RED, COL_FINAL], alpha=0.85, width=0.5)
    for b, v in zip(bars, nhigh):
        ax2.annotate(f"{v}", (b.get_x()+b.get_width()/2, v+0.4), ha="center",
                     fontsize=11, fontweight="bold")
    ax2.set_ylabel("APR-Score high-risk windows (nHigh, of 69)")
    ax2.set_ylim(0, 38)
    ax2.set_title("(b) De-alanine redesign reduces aggregation windows")
    ax2.text(0.02, 0.95, "A% 38.8 -> 36.2; ipTM 0.84 -> 0.937\ndG -89.7 -> -129.9 kcal/mol",
             transform=ax2.transAxes, fontsize=8.5, color=COL_ACC, va="top")
    fig.suptitle("Quantified effects of the developability-integrated design loop",
                 fontsize=11.5)
    save(fig, "Fig5_developability")

# =====================================================================
# Figure 6. Dimer modeling: V4 bridges the VEGF dimer interface
# =====================================================================
def fig6_dimer():
    cands = ["V1", "V2", "V4"]
    # std-sequence recheck (2026-08-18)
    ab = np.array([0.610, 0.581, 0.619])
    ac = np.array([0.117, 0.075, 0.682])
    bc = np.array([0.120, 0.073, 0.662])
    x = np.arange(3); w = 0.26

    fig, ax = plt.subplots(figsize=(7.4, 4.0))
    ax.bar(x - w, ab, w, label="VEGF-A to VEGF-B (A,B)", color="#9db8d2")
    ax.bar(x, ac, w, label="VEGF-A to binder (A,C)", color=COL_FINAL)
    ax.bar(x + w, bc, w, label="VEGF-B to binder (B,C)", color=COL_ACC)
    ax.set_xticks(x); ax.set_xticklabels(cands)
    ax.set_ylabel("OpenFold3 chain-pair ipTM")
    ax.set_ylim(0, 1.0)
    ax.axhline(0.5, color=COL_RED, linestyle="--", linewidth=0.8)
    ax.text(1.0, 0.53, "0.5 (binder–VEGF indicative)", fontsize=8, color=COL_RED, ha="center", va="bottom")
    for xi in range(3):
        ax.annotate(f"{ac[xi]:.2f}", (xi, ac[xi]+0.02), ha="center", fontsize=8.5,
                    color=COL_FINAL, fontweight="bold")
    ax.set_title("VEGF dimer (2xVEGF165) + binder template modeling: V4 bridges the dimer interface\n(single best-seed OF3 run per candidate)")
    ax.legend(frameon=False, fontsize=8)
    ax.text(0.02, 0.95, "V4 scores ~0.68 against both VEGF chains\n-> potential dimerization-blocking mechanism\n(requires SPR/BLI confirmation)",
            transform=ax.transAxes, fontsize=8.5, color=COL_ACC, va="top")
    save(fig, "Fig6_dimer")

if __name__ == "__main__":
    fig1_pipeline()
    fig2_redesign_dG()
    fig3_full_length()
    fig4_crossval()
    fig5_developability()
    fig6_dimer()
    print("ALL FIGURES DONE")
