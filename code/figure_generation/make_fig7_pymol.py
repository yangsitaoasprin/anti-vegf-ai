# -*- coding: utf-8 -*-
# =============================================================================
# Author : Yang Sitao (杨四涛)
# Affiliation : School of Pharmaceutical Sciences, Dali University;
#               The Third Affiliated Hospital of Dali University
# Contact : yangsitaoasprin@swpu.edu.cn
# License : MIT (see ../../LICENSE)
# Purpose: Reproduces Fig. 2 of the manuscript (PyMOL rendering of the three
#          NOTE: this file keeps its legacy name; the manuscript's figure
#          NUMBER is 2 (build_docx.py maps "2" -> Fig7_binding_pose.png)
#          binder-VEGF-A RBD complexes).
# =============================================================================
"""
Fig. 2 (file name Fig7_binding_pose) — publication-grade rendering with PyMOL
(open-source, headless ray tracing).

For each of the three redesigned binders (V1/V2/V4 Boltz-1 RBD complexes):
  - VEGF-A RBD (chain A) drawn as grey cartoon
  - miniprotein binder (chain B) drawn as coloured cartoon
  - interface residues (< 5 A heavy-atom contact across chains) shown as
    red (VEGF side) / orange (binder side) sticks + spheres so the binding
    patch pops out
  - depth-cued lighting + ray tracing (anti-aliased, AO-ish) -> journal quality
Each panel is ray-traced to a high-res PNG, then montaged into one figure.

Run with the pymol env:
  conda run -n pymol python make_fig7_pymol.py <OUT_DIR>
"""
import json
import os, sys, glob
import numpy as np
from PIL import Image, ImageDraw, ImageFont

try:
    from pymol import cmd as _PYMOL_CMD
except Exception:
    _PYMOL_CMD = None

def _safe_set(name, val, sel="all"):
    """Set a PyMOL setting, silently ignoring unknown-setting errors that vary
    across PyMOL builds (e.g. open-source PyMOL lacks 'fog_end')."""
    cmd = _PYMOL_CMD
    if cmd is None:
        try:
            from pymol import cmd as c
        except Exception:
            return
        cmd = c
    try:
        cmd.set(name, val, sel)
    except Exception:
        pass

# Representative structures ship inside this repository (repro/data/representative_structures/).
# Locate the repo root by walking up from this script (the data dir lives at
# repro/data/, i.e. three levels above repro/code/figure_generation/).
_HERE = os.path.dirname(os.path.abspath(__file__))


def _find_repo_root(start):
    d = start
    for _ in range(6):
        if os.path.isdir(os.path.join(d, "data", "representative_structures")):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    return os.path.dirname(os.path.dirname(start))


_REPO_ROOT = _find_repo_root(_HERE)                # .../repro
REP = os.path.join(_REPO_ROOT, "data", "representative_structures")
OUT  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.getcwd(), "fig7_out")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------------------
# Interface-residue counts (printed under each panel)
# ---------------------------------------------------------------------------
# The Fig. 2 caption promises "numbers below each panel"; these must come from
# the SAME structure that is drawn and use the SAME definition as Fig. 3, so a
# later edit cannot let the two figures drift apart.
IFACE_CUTOFF = 5.0


def parse_pdb(path):
    """{chain: {resnum: [(x, y, z), ...]}} for heavy (non-H) atoms."""
    chains = {}
    with open(path) as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if line[12:16].strip().startswith("H"):
                continue
            try:
                resnum = int(line[22:26])
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            chains.setdefault(line[21], {}).setdefault(resnum, []).append((x, y, z))
    return chains


def interface_counts(chains, cutoff=IFACE_CUTOFF):
    """(vegf_side, binder_side, n_vegf_res, n_binder_res).

    A residue is an interface residue when ANY of its heavy atoms is within
    `cutoff` of ANY heavy atom of the partner chain.  NEVER break out of the
    inner loop early: doing so makes each chain-A atom claim at most one
    chain-B residue and under-counts the binder side (V1 16 instead of 25,
    V2 25 instead of 29, V4 14 instead of 17 -- the bug fixed on 2026-09-15
    in make_fig8_interface_stats.py)."""
    A = chains.get("A", {})   # VEGF-A RBD
    B = chains.get("B", {})   # miniprotein binder
    cut2 = cutoff * cutoff
    vi, bi = set(), set()
    for ra, aa in A.items():
        for rb, ab in B.items():
            if any((xa - xb) ** 2 + (ya - yb) ** 2 + (za - zb) ** 2 < cut2
                   for xa, ya, za in aa for xb, yb, zb in ab):
                vi.add(ra); bi.add(rb)
    return len(vi), len(bi), len(A), len(B)


# (candidate key, label, pdb path, binder colour hex without '#')
CANDS = [
    ("V1", "V1  #1_T0.3_s132",
     os.path.join(REP, "V1_s132_model_0.pdb"),
     "4f8bc9"),   # blue
    ("V2", "V2  #2c_T0.1_s18",
     os.path.join(REP, "V2_s18_model_0.pdb"),
     "5a9e6f"),   # green
    ("V4", "V4  T0.1_s35",
     os.path.join(REP, "V4_s35_model_0.pdb"),
     "e08a3c"),   # orange
]

def render_panel(pdb, binder_hex, out_png, azim=140, elev=18, d=0):
    import pymol
    from pymol import cmd
    pymol.finish_launching(["pymol", "-cq"])  # quiet headless
    cmd.reinitialize()
    cmd.load(pdb, "C")
    # Chain convention (matches the manuscript caption and ALL representative PDBs):
    # VEGF-A RBD is on chain A, the miniprotein binder on chain B. NOTE: the binder can
    # be LONGER than the 98-aa RBD (e.g. V2 is 100 aa), so a length-based heuristic would
    # invert V2's colours (grey VEGF / green binder). We therefore assign by chain ID.
    chains = cmd.get_chains("C")
    if "A" in chains and "B" in chains:
        veg_ch, bin_ch = "A", "B"
    else:
        lens = {ch: cmd.count_atoms(f"C and chain {ch}") for ch in chains}
        veg_ch = max(lens, key=lens.get)
        bin_ch = min(lens, key=lens.get)

    cmd.hide("everything")
    # ----- VEGF cartoon (grey) -----
    cmd.show("cartoon", f"C and chain {veg_ch}")
    cmd.color("grey80", f"C and chain {veg_ch}")
    # ----- binder cartoon (colour) -----
    cmd.show("cartoon", f"C and chain {bin_ch}")
    cmd.color(f"0x{binder_hex}", f"C and chain {bin_ch}")
    # make cartoon slightly fat + smooth
    cmd.set("cartoon_flat_sheets", 0)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("cartoon_loop_radius", 0.4)

    # ----- interface residues across chains (< 5 A) -----
    cmd.select("ifA", f"(chain {veg_ch} within 5 of chain {bin_ch})")
    cmd.select("ifB", f"(chain {bin_ch} within 5 of chain {veg_ch})")
    # receptor-side interface: red sticks + spheres
    cmd.show("sticks", "ifA")
    cmd.show("spheres", "ifA")
    cmd.color("red", "ifA")
    cmd.set("sphere_scale", 0.45, "ifA")
    cmd.set("stick_radius", 0.18, "ifA")
    # binder-side interface: orange sticks + spheres
    cmd.show("sticks", "ifB")
    cmd.show("spheres", "ifB")
    cmd.color("orange", "ifB")
    cmd.set("sphere_scale", 0.45, "ifB")
    cmd.set("stick_radius", 0.18, "ifB")
    # dim the non-interface sticks of the binder a touch by valence colour? leave cartoon colour
    cmd.show("cartoon", "ifA")  # keep cartoon underneath
    cmd.set("cartoon_color", "grey80", "ifA")
    cmd.set("cartoon_color", f"0x{binder_hex}", "ifB")

    # ----- lighting / background -----
    cmd.bg_color("white")
    _safe_set("depth_cue", 1)
    _safe_set("fog_start", 0.35)
    _safe_set("fog_end", 0.95)
    _safe_set("light_count", 3)
    _safe_set("specular", 0.45)
    _safe_set("shininess", 40)
    _safe_set("ambient", 0.30)
    _safe_set("direct", 0.70)
    _safe_set("ray_shadows", 1)
    _safe_set("antialias", 2)
    _safe_set("ray_trace_mode", 1)   # normal-blended shading for crisp protein

    # ----- orient & zoom on interface -----
    cmd.select("iface_union", "ifA or ifB")
    cmd.orient("iface_union")
    cmd.zoom("iface_union", 18)     # padding around the binding patch
    cmd.turn("y", azim)
    cmd.turn("x", elev)
    if d:
        cmd.move("z", d)

    cmd.ray(1400, 1400)
    cmd.png(out_png, 1400, 1400, dpi=300, ray=1)
    cmd.reinitialize()
    return out_png

# ---------------------------------------------------------------------------
# Figure geometry and typography.
#
# paper_check.py couples the shipped PNG back to these constants.
#
# Typography is expressed in POINTS, never in pixels.  Elsevier requires the
# lettering inside artwork to end up at >= 7 pt at the printed size -- "Artwork
# sizing" (https://www.elsevier.com/about/policies-and-standards/author/
# artwork-and-media-instructions/artwork-sizing): "the lettering on the artwork
# should have a finished, printed size of 7 pt for normal text and no smaller
# than 6 pt for subscript and superscript characters. Smaller lettering will
# yield text that is hardly legible."  The manuscript places this figure at
# FIG_WIDTH_IN inches (build_docx.FIG_WIDTH["2"]), so the canvas carries
# PANEL_W*3/(FIG_WIDTH_IN*72) pixels per point.  A pixel-hardcoded size thus
# silently rescales whenever either number changes: the previous 46/32/28 px
# annotation printed at 5.1/3.6/3.1 pt and the legend was unreadable.
# ---------------------------------------------------------------------------
PANEL_W = 1400                       # ray-traced panel edge (render_panel)
FIG_WIDTH_IN = 6.5                   # must equal build_docx.FIG_WIDTH["2"]
PX_PER_PT = (PANEL_W * 3.0) / (FIG_WIDTH_IN * 72.0)      # ~8.97 px per point
FS_PT = 7.0                          # every annotation in this figure
GAP_LABEL = 30                       # px of white around the panel label
GAP_COUNTS = 30                      # px of white around the counts line
GAP_LEGEND = 46                      # px of white below/above the legend row
                                     # (kept larger than GAP_COUNTS so the
                                     # full-width legend reads as a footer and
                                     # not as an annotation of the nearest panel)
DOT_TEXT_GAP = 10                    # px between a legend dot and its label
LEGEND_ITEM_FRAC = 2.5               # inter-item gap, in label-font sizes
EDGE_MARGIN = 40                     # px kept clear at the canvas left/right

# Legacy geometry of the shipped 2026-09-15 13:29 render.  Only the lossless
# slicer (论文/_fig2_remontage2.py) needs these, to crop the pre-counts PNG.
PAD_TOP_LEGACY = 90
PAD_BOT_LEGACY = 70


def _pt(v):
    """Points -> pixels on the shipped canvas."""
    return int(round(v * PX_PER_PT))


def _font_path(font):
    """Filesystem path of the resolved face, or "" when it is not a real file.

    ImageFont.load_default() returns a face backed by an in-memory BytesIO, so
    this returns "" for the silent bitmap fallback -- which is how the gate
    detects that the fallback happened at all.  It also keeps the layout JSON
    serialisable, since a BytesIO is not."""
    p = getattr(font, "path", None)
    return p if isinstance(p, str) else ""


def _text_h(d, txt, font):
    """(ink height, ink width) of `txt` in `font`; textbbox is offset so it is
    [top, bottom] that matter, not the nominal font size."""
    b = d.textbbox((0, 0), txt, font=font)
    return b[3] - b[1], b[2] - b[0]


def _font(size, bold=False):
    """Resolve a scalable font.

    ImageFont.truetype() RAISES on a bare name when the face is not installed,
    and the old code caught that and fell back to the ~11-px bitmap default --
    which silently shrank the panel labels from 46 px to 9 px, five times
    smaller than the published figure, while leaving a perfectly plausible
    looking image behind.  matplotlib bundles DejaVu (both the house style and
    present on Windows and in the WSL pymol env), so try that first, then the
    platform faces, and only then a scalable default."""
    names = ["DejaVuSans-Bold.ttf", "DejaVuSans.ttf"] if bold else ["DejaVuSans.ttf"]
    paths = []
    try:
        import matplotlib
        _mpl = os.path.join(os.path.dirname(matplotlib.__file__),
                            "mpl-data", "fonts", "ttf")
        paths += [os.path.join(_mpl, n) for n in names]
    except Exception:
        pass
    paths += list(names)
    paths += ["/usr/share/fonts/truetype/dejavu/" + n for n in names]
    paths += ["C:/Windows/Fonts/arialbd.ttf", "C:/Windows/Fonts/arial.ttf"]
    for p in paths:
        try:
            return ImageFont.truetype(p, size)
        except Exception:
            continue
    try:
        return ImageFont.load_default(size)      # Pillow >= 10.1 is scalable
    except Exception:
        return ImageFont.load_default()


LEGEND = [("VEGF-A RBD", "cccccc"), ("miniprotein binder", "4f8bc9"),
          ("interface (VEGF side)", "ff0000"), ("interface (binder side)", "ffa500")]
# widest plausible counts line -- used only to size the band, never drawn
COUNTS_TMPL = "interface residues: VEGF 88 / binder 88"


def montage(panels, labels, out_png, counts=None):
    """Combine 3 square PNGs into one labelled figure with a legend.

    `counts` is an optional list of (vegf_side, binder_side) interface-residue
    counts, one per panel.  When given, a counts line is drawn under each panel
    -- the manuscript Fig. 2 caption promises exactly that, so the figure must
    not be assembled without them.

    Every band is derived from the measured text metrics, never from a
    hardcoded strip height, so changing FS_PT can neither clip an annotation nor
    leave a stale band of white space behind.  The realised layout is written to
    <out_png>_layout.json so the gate can verify it after the fact."""
    W = PANEL_W
    flab = _font(_pt(FS_PT), bold=True)
    fcnt = _font(_pt(FS_PT))
    fleg = _font(_pt(FS_PT))
    probe = ImageDraw.Draw(Image.new("RGB", (8, 8)))

    lab_h = max(_text_h(probe, l, flab)[0] for l in labels)
    cnt_h = _text_h(probe, COUNTS_TMPL, fcnt)[0] if counts else 0
    leg_bb = probe.textbbox((0, 0), "Ag", font=fleg)     # ascender + descender
    dot_r = max(4, int(round(_pt(FS_PT) / 2.0)))

    pad_top = lab_h + 2 * GAP_LABEL
    pad_counts = (cnt_h + 2 * GAP_COUNTS) if counts else 0
    pad_bot = (leg_bb[3] - leg_bb[1]) + 2 * GAP_LEGEND
    H = pad_top + W + pad_counts + pad_bot

    canvas = Image.new("RGB", (W * 3, H), "white")
    d = ImageDraw.Draw(canvas)
    ink = {}

    for i, (p, lab) in enumerate(zip(panels, labels)):
        im = Image.open(p).convert("RGB")
        im = im.resize((W, W), Image.LANCZOS)
        canvas.paste(im, (i * W, pad_top))
        # label, centred over its own panel
        lb = d.textbbox((0, 0), lab, font=flab)
        ly = pad_top - GAP_LABEL - lb[3]
        d.text(((i * W) + (W - (lb[2] - lb[0])) / 2.0, ly), lab, font=flab,
               fill="black")
        if i == 0:
            ink["label"] = [ly + lb[1], ly + lb[3]]
        # interface-residue counts -- the caption's promise, same centring
        if counts:
            vs, bs = counts[i]
            ctxt = "interface residues: VEGF %d / binder %d" % (vs, bs)
            cb = d.textbbox((0, 0), ctxt, font=fcnt)
            cy = pad_top + W + GAP_COUNTS - cb[1]
            d.text(((i * W) + (W - (cb[2] - cb[0])) / 2.0, cy), ctxt, font=fcnt,
                   fill="black")
            if i == 0:
                ink["counts"] = [cy + cb[1], cy + cb[3]]

    # legend, centred across the whole canvas
    item_gap = int(round(LEGEND_ITEM_FRAC * _pt(FS_PT)))
    tw = [d.textbbox((0, 0), t, font=fleg)[2] for t, _c in LEGEND]
    while True:
        total = (sum(tw) + len(LEGEND) * (2 * dot_r + DOT_TEXT_GAP)
                 + (len(LEGEND) - 1) * item_gap)
        if total <= W * 3 - 2 * EDGE_MARGIN or item_gap <= 16:
            break
        item_gap -= 4
    x = (W * 3 - total) / 2.0
    ly = H - GAP_LEGEND - leg_bb[3]
    cy = ly + (leg_bb[1] + leg_bb[3]) / 2.0
    for (txt, hexc), w in zip(LEGEND, tw):
        d.ellipse([x, cy - dot_r, x + 2 * dot_r, cy + dot_r], fill="#" + hexc)
        d.text((x + 2 * dot_r + DOT_TEXT_GAP, ly), txt, font=fleg, fill="black")
        x += 2 * dot_r + DOT_TEXT_GAP + w + item_gap
    ink["legend"] = [ly + leg_bb[1], ly + leg_bb[3]]

    canvas.save(out_png)

    with open(os.path.splitext(out_png)[0] + "_layout.json", "w") as f:
        json.dump({
            "fig_width_in": FIG_WIDTH_IN, "panel_w": W,
            "px_per_pt": round(PX_PER_PT, 4), "fs_pt": FS_PT,
            "pad_top": pad_top, "pad_counts": pad_counts, "pad_bot": pad_bot,
            "size": [W * 3, H],
            "bands": {"panel": [pad_top, pad_top + W],
                      "counts": [pad_top + W, pad_top + W + pad_counts],
                      "legend": [H - pad_bot, H]},
            "ink": ink,
            "fonts": {k: _font_path(f)
                      for k, f in [("label", flab), ("counts", fcnt),
                                   ("legend", fleg)]},
        }, f, indent=2)
    return out_png

def main():
    panels, labels, counts = [], [], []
    for key, lab, pdb, hexc in CANDS:
        op = os.path.join(OUT, f"_panel_{key}.png")
        # slightly different azimuths so the 3 panels don't look identical
        az = {0: 140, 1: 120, 2: 160}[len(panels)]
        render_panel(pdb, hexc, op, azim=az, elev=18)
        panels.append(op)
        labels.append(lab)
        # counts come from the SAME structure that was just drawn, with the
        # same 5 A heavy-atom definition as Fig. 3 -- never hardcoded
        vs, bs, _na, _nb = interface_counts(parse_pdb(pdb))
        counts.append((vs, bs))
    # write a manifest so the montage step (run under base python) can find panels
    with open(os.path.join(OUT, "_panels_manifest.txt"), "w") as f:
        for p, lab in zip(panels, labels):
            f.write(f"{p}\t{lab}\n")
    # assemble the 3 panels into the final figure (Fig. 2 in the manuscript)
    final = os.path.join(OUT, "Fig7_binding_pose.png")
    montage(panels, labels, final, counts=counts)
    with open(os.path.join(OUT, "Fig7_binding_pose_iface_counts.json"), "w") as f:
        json.dump({
            "definition": "heavy-atom contact < %.1f A; chain A = VEGF-A RBD,"
                          " chain B = binder" % IFACE_CUTOFF,
            "rows": {k: {"vegf_side": v, "binder_side": b}
                     for (k, _l, _p, _h), (v, b) in zip(CANDS, counts)},
        }, f, indent=2)
    print("panels rendered:", panels)
    print("interface counts:", counts)
    print("figure written:", final)

if __name__ == "__main__":
    main()
