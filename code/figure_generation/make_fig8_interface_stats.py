# -*- coding: utf-8 -*-
"""
make_fig8_interface_stats.py
============================
Generator for Manuscript Fig. 3 / figures/Fig8_interface_stats
("Interface statistics of the three redesigned binders on the VEGF-A RBD").

This script REPLACES the previously hand-built static SVG (which had the
interface counts hard-coded and therefore was not reproducible). It now:

  1. Parses the three final-candidate Boltz-1 representative complexes
     (chain A = VEGF-A RBD, 98 residues; chain B = redesigned binder).
  2. Counts interface residues on each side as residues that have at least
     one heavy atom within 5.0 A of any heavy atom of the partner chain
     ("< 5 A heavy-atom contact", per the manuscript caption).
  3. Computes the binder-interface fraction = binder_interface / binder_length.
  4. Renders the two-panel figure (a) interface size, (b) binder fraction
     and writes PNG + SVG alongside the other manuscript figures.

Run:  python make_fig8_interface_stats.py [FIG_DIR]
If FIG_DIR is omitted the script defaults to the folder containing this file
(../论文/figures).  Numbers are also written to Fig8_interface_stats_numbers.json
for traceability / verification against the manuscript text.
"""

import os
import sys
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Printable-minimum text guard: matplotlib sizes text relative to the FIGURE, so
# a figure authored wider than it is placed prints its labels too small (this
# one is authored 9.0 in wide, placed 6.0 in).  See figstyle.py; it calibrates
# against FIG_WIDTH in build_docx.py, the single source of truth.
import figstyle

# ---------------------------------------------------------------------------
# Input: final-candidate representative Boltz-1 complexes (truncated VEGF RBD)
# ---------------------------------------------------------------------------
# Resolved by walking up from this file rather than hard-coded, so the same
# script works both in the authoring tree (论文/make_fig8_interface_stats.py,
# data under <root>/repro/data/...) and inside the published repro/ bundle
# (repro/code/figure_generation/make_fig8_interface_stats.py, data under
# <repro>/data/...).
def _find_rep_dir():
    here = os.path.dirname(os.path.abspath(__file__))
    for _ in range(6):
        for parts in (("repro", "data", "representative_structures"),
                      ("data", "representative_structures")):
            cand = os.path.join(here, *parts)
            if os.path.isdir(cand):
                return cand
        parent = os.path.dirname(here)
        if parent == here:
            break
        here = parent
    return None


REP_DIR = _find_rep_dir()
if REP_DIR is None:
    raise SystemExit(
        "make_fig8_interface_stats.py: could not locate "
        "'data/representative_structures' by walking up from %s"
        % os.path.dirname(os.path.abspath(__file__)))
print("input dir:", REP_DIR, flush=True)

PDBS = {
    "V1": os.path.join(REP_DIR, "V1_s132_model_0.pdb"),   # #1_T0.3_s132  (80 aa binder)
    "V2": os.path.join(REP_DIR, "V2_s18_model_0.pdb"),    # #2c_T0.1_s18 (100 aa binder)
    "V4": os.path.join(REP_DIR, "V4_s35_model_0.pdb"),    # T0.1_s35     (80 aa binder)
}

CONTACT_CUTOFF = 5.0  # Angstrom, heavy-atom, per manuscript Fig.3 caption

# Colours copied from the original static SVG so the new figure is visually
# faithful to the previously published layout.
C_VEGF = "#d43f3f"   # VEGF-side bars (panel a)
C_BIND = "#e8a33d"   # binder-side bars (panel a)
C_FRAC = {"V1": "#2a6fbd", "V2": "#2e8b57", "V4": "#e07b39"}  # panel b


# ---------------------------------------------------------------------------
# PDB parsing + interface counting
# ---------------------------------------------------------------------------
def parse_pdb(path):
    """Return {chain: {resnum: [(x,y,z), ...]}} for heavy (non-H) atoms."""
    chains = {}
    with open(path) as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom = line[12:16].strip()
            if atom.startswith("H"):
                continue
            chain = line[21]
            resnum = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            chains.setdefault(chain, {}).setdefault(resnum, []).append((x, y, z))
    return chains


def interface_counts(chains, cutoff=CONTACT_CUTOFF):
    """Count interface residues on each side (heavy-atom contact < cutoff A).

    A residue counts as an interface residue when ANY of its heavy atoms lies
    within `cutoff` of ANY heavy atom of the partner chain, and EVERY partner
    residue satisfying this is collected.

    BUGFIX 2026-09-15: the previous implementation broke out of the inner loop
    as soon as it found the FIRST partner residue contacting a given A atom
    (`bind_iface.add(rb); break`).  Each A atom therefore "claimed" at most one
    B residue, so the binder side was systematically under-counted
    (V1 16 instead of 25, V2 25 instead of 29, V4 14 instead of 17) while the
    VEGF side - which only needs a boolean - was unaffected.  Those wrong
    numbers were what put the recomputed figure into figures/quarantine/.
    The correct counts reproduce the published Fig. 3 (18/25, 27/29, 19/17).
    """
    A = chains["A"]  # VEGF-A RBD
    B = chains["B"]  # binder
    cut2 = cutoff * cutoff
    vegf_iface = set()
    bind_iface = set()
    for ra, atoms_a in A.items():
        for rb, atoms_b in B.items():
            if any((xa - xb) ** 2 + (ya - yb) ** 2 + (za - zb) ** 2 < cut2
                   for xa, ya, za in atoms_a
                   for xb, yb, zb in atoms_b):
                # this residue pair is in contact: record BOTH sides and keep
                # scanning - never break out early
                vegf_iface.add(ra)
                bind_iface.add(rb)
    return len(vegf_iface), len(bind_iface), len(A), len(B)


def compute_all():
    rows = {}
    for name, pdb in PDBS.items():
        if not os.path.exists(pdb):
            raise FileNotFoundError(f"Missing representative PDB for {name}: {pdb}")
        chains = parse_pdb(pdb)
        vegf_n, bind_n, vegf_len, bind_len = interface_counts(chains)
        rows[name] = {
            "vegf_side": vegf_n,
            "binder_side": bind_n,
            "vegf_len": vegf_len,
            "binder_len": bind_len,
            "binder_frac": round(bind_n / bind_len * 100.0, 1),
        }
    return rows


# ---------------------------------------------------------------------------
# Figure rendering (mirrors the original static SVG layout)
# ---------------------------------------------------------------------------
def render(rows, fig_dir):
    order = ["V1", "V2", "V4"]
    vegf = [rows[k]["vegf_side"] for k in order]
    bind = [rows[k]["binder_side"] for k in order]
    frac = [rows[k]["binder_frac"] for k in order]

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(9.0, 3.8), gridspec_kw={"width_ratios": [1, 1]}
    )
    fig.subplots_adjust(left=0.07, right=0.98, wspace=0.32, bottom=0.22, top=0.82)

    # ---- panel (a): interface size -----------------------------------------
    x = range(len(order))
    w = 0.38
    ax1.bar([i - w / 2 for i in x], vegf, width=w, color=C_VEGF, alpha=0.85,
            label="VEGF-side")
    ax1.bar([i + w / 2 for i in x], bind, width=w, color=C_BIND, alpha=0.9,
            label="binder-side")
    for i, (v, b) in enumerate(zip(vegf, bind)):
        ax1.text(i - w / 2, v + 0.6, str(v), ha="center", va="bottom",
                 fontsize=9)
        ax1.text(i + w / 2, b + 0.6, str(b), ha="center", va="bottom",
                 fontsize=9)
    ax1.set_xticks(list(x))
    ax1.set_xticklabels(order, fontsize=10)
    ax1.set_ylabel("interface residues (< 5 Å)", fontsize=9)
    ax1.set_ylim(0, max(vegf + bind) * 1.18)
    ax1.set_title("(a) Interface size", fontsize=10, loc="left")
    ax1.legend(handles=[Patch(facecolor=C_VEGF, alpha=0.85, label="VEGF-side"),
                        Patch(facecolor=C_BIND, alpha=0.9, label="binder-side")],
               loc="upper left", fontsize=8, frameon=False)

    # ---- panel (b): binder-interface fraction ------------------------------
    bars = ax2.bar(list(x), frac, width=0.55,
                   color=[C_FRAC[k] for k in order], alpha=0.9)
    for i, f in enumerate(frac):
        ax2.text(i, f + 0.8, f"{f:.0f}%", ha="center", va="bottom", fontsize=9)
    ax2.set_xticks(list(x))
    ax2.set_xticklabels(order, fontsize=10)
    ax2.set_ylabel("binder interface residues (% of binder)", fontsize=9)
    ax2.set_ylim(0, max(frac) * 1.25)
    ax2.set_title("(b) Binder-interface fraction", fontsize=10, loc="left")

    fig.suptitle(
        "Interface statistics of the three redesigned binders on the VEGF-A RBD "
        "(Boltz-1 complexes)",
        fontsize=9.5, y=0.96,
    )

    png = os.path.join(fig_dir, "Fig8_interface_stats.png")
    svg = os.path.join(fig_dir, "Fig8_interface_stats.svg")
    # Enforce the 7 pt printed-size floor (see figstyle.py) before writing.
    rep = figstyle.normalize(fig, "Fig8_interface_stats.png", out_dir=fig_dir)
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    print("saved Fig8_interface_stats |", figstyle.report_line(rep), flush=True)
    for ov in rep["overlaps_after"]:
        print("   !! text overlap %.2f: %r <> %r" % (ov[3], ov[0], ov[1]),
              flush=True)
    return png, svg


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    fig_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    rows = compute_all()
    png, svg = render(rows, fig_dir)

    numbers = {
        "definition": f"heavy-atom contact < {CONTACT_CUTOFF} A, per-chain interface residue counts",
        "rows": rows,
    }
    json_path = os.path.join(fig_dir, "Fig8_interface_stats_numbers.json")
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(numbers, fh, indent=2)

    # human-readable summary to stdout-equivalent log
    out = os.path.join(fig_dir, "Fig8_interface_stats_run.log")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("Fig8 interface statistics (computed from representative PDBs)\n")
        fh.write(f"definition: {numbers['definition']}\n\n")
        fh.write(f"{'cand':<5}{'VEGF-side':>10}{'binder-side':>13}{'binder_len':>12}{'frac_%':>8}\n")
        for k in ["V1", "V2", "V4"]:
            r = rows[k]
            fh.write(f"{k:<5}{r['vegf_side']:>10}{r['binder_side']:>13}"
                     f"{r['binder_len']:>12}{r['binder_frac']:>8.1f}\n")
        fh.write(f"\nwrote: {png}\nwrote: {svg}\nwrote: {json_path}\n")

    print(f"OK wrote {png}\nOK wrote {svg}\nOK wrote {json_path}")
    for k in ["V1", "V2", "V4"]:
        r = rows[k]
        print(f"{k}: VEGF-side={r['vegf_side']} binder-side={r['binder_side']} "
              f"frac={r['binder_frac']}%")


if __name__ == "__main__":
    main()
