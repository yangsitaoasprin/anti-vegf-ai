# -*- coding: utf-8 -*-
# =============================================================================
# Author : Yang Sitao (杨四涛)
# Affiliation : School of Pharmaceutical Sciences, Dali University;
#               The Third Affiliated Hospital of Dali University
# Contact : yangsitaoasprin@126.com
# License : MIT (see ../../LICENSE)
# Purpose: Reproduces Fig. 7 of the manuscript (PyMOL rendering of the three
#          binder-VEGF-A RBD complexes).
# =============================================================================
"""
Fig.7 — publication-grade rendering with PyMOL (open-source, headless ray tracing).

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
import os, sys, glob
import numpy as np
from PIL import Image, ImageDraw, ImageFont

# Representative structures ship inside this repository (data/representative_structures/).
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # .../repro
REP = os.path.join(_REPO_ROOT, "data", "representative_structures")
OUT  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.getcwd(), "fig7_out")
os.makedirs(OUT, exist_ok=True)

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
    # robust chain detection: longer chain = VEGF RBD (chain A), shorter = binder (chain B)
    chains = cmd.get_chains("C")
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
    cmd.set("depth_cue", 1)
    cmd.set("fog_start", 0.35)
    cmd.set("fog_end", 0.95)
    cmd.set("light_count", 3)
    cmd.set("specular", 0.45)
    cmd.set("shininess", 40)
    cmd.set("ambient", 0.30)
    cmd.set("direct", 0.70)
    cmd.set("ray_shadows", 1)
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 1)   # normal-blended shading for crisp protein

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

def montage(panels, labels, out_png):
    """Combine 3 square PNGs into one labelled figure with a legend."""
    W = 1400
    pad_top = 90
    pad_bot = 70
    H = W + pad_top + pad_bot
    canvas = Image.new("RGB", (W * 3, H), "white")
    d = ImageDraw.Draw(canvas)
    try:
        fbig  = ImageFont.truetype("DejaVuSans-Bold.ttf", 46)
        fsmall= ImageFont.truetype("DejaVuSans.ttf", 30)
    except Exception:
        fbig  = ImageFont.load_default()
        fsmall= ImageFont.load_default()
    for i, (p, lab) in enumerate(zip(panels, labels)):
        im = Image.open(p).convert("RGB")
        im = im.resize((W, W), Image.LANCZOS)
        canvas.paste(im, (i * W, pad_top))
        # label
        tb = d.textbbox((0, 0), lab, font=fbig)
        tw = tb[2] - tb[0]
        d.text(((i * W) + (W - tw) / 2, pad_top - 66), lab, font=fbig, fill="black")
    # legend at bottom
    leg = [("VEGF-A RBD", "grey80"), ("miniprotein binder", "4f8bc9"),
           ("interface (VEGF side)", "red"), ("interface (binder side)", "orange")]
    x = 60
    y = H - 42
    try:
        fleg = ImageFont.truetype("DejaVuSans.ttf", 28)
    except Exception:
        fleg = ImageFont.load_default()
    for txt, hexc in leg:
        d.ellipse([x, y - 14, x + 28, y + 14], fill=("#" + hexc))
        d.text((x + 38, y - 12), txt, font=fleg, fill="black")
        x += 38 + d.textbbox((0, 0), txt, font=fleg)[2] + 70
    canvas.save(out_png)
    return out_png

def main():
    panels, labels = [], []
    for key, lab, pdb, hexc in CANDS:
        op = os.path.join(OUT, f"_panel_{key}.png")
        # slightly different azimuths so the 3 panels don't look identical
        az = {0: 140, 1: 120, 2: 160}[len(panels)]
        render_panel(pdb, hexc, op, azim=az, elev=18)
        panels.append(op)
        labels.append(lab)
    # write a manifest so the montage step (run under base python) can find panels
    with open(os.path.join(OUT, "_panels_manifest.txt"), "w") as f:
        for p, lab in zip(panels, labels):
            f.write(f"{p}\t{lab}\n")
    print("panels rendered:", panels)

if __name__ == "__main__":
    main()
