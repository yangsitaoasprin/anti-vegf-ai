# -*- coding: utf-8 -*-
# =============================================================================
# Author : Yang Sitao (杨四涛)
# Affiliation : School of Pharmaceutical Sciences, Dali University;
#               The Third Affiliated Hospital of Dali University
# Contact : yangsitaoasprin@swpu.edu.cn
# License : MIT (see ../../LICENSE)
# Purpose: Reproduces the two PyMOL panels embedded in the Graphical Abstract.
# =============================================================================
"""
Render the two outcome panels of the Graphical Abstract with PyMOL (headless,
ray-traced), reusing the Fig.7 visual style (grey antigen cartoon, coloured
binder cartoon, red/orange interface sticks + spheres, white bg, AO-ish shading).

  GA_V1_contact.png : V1 weakly contacts VEGF-A chain A (Boltz monomer complex)
  GA_V4_bridge.png  : V4 bridges the VEGF-A dimer (chain A + B) (OF3 dimer complex)

Run under the WSL pymol conda env:
  /root/miniconda3/envs/pymol/bin/python render_ga_outcome.py
"""
import os
import pymol
from pymol import cmd

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
OUT = _HERE                                        # write PNGs next to this script
os.makedirs(OUT, exist_ok=True)

V1_PDB = os.path.join(REP, "V1_s132_model_0.pdb")
V4_CIF = os.path.join(REP, "V4_dimer_bridge_model.cif")


def render_complex(pdb, binder_hex, out_png, azim=140, elev=18):
    """Render one binder+VEGF complex. Binder = shortest chain by atom count;
    all other chains are treated as the antigen (VEGF, grey)."""
    pymol.finish_launching(["pymol", "-cq"])
    cmd.reinitialize()
    obj = "C"
    cmd.load(pdb, obj)

    chains = cmd.get_chains(obj)
    lens = {ch: cmd.count_atoms(f"{obj} and chain {ch}") for ch in chains}
    print("  chains:", lens)
    bin_ch = min(lens, key=lens.get)              # shortest chain = miniprotein binder
    veg_chs = [ch for ch in chains if ch != bin_ch]
    print("  binder =", bin_ch, "| antigen chains =", veg_chs)

    cmd.hide("everything")
    veg_sel = " or ".join(f"chain {ch}" for ch in veg_chs)
    # antigen cartoon (grey)
    cmd.show("cartoon", f"{obj} and ({veg_sel})")
    cmd.color("grey80", f"{obj} and ({veg_sel})")
    # binder cartoon (colour)
    cmd.show("cartoon", f"{obj} and chain {bin_ch}")
    cmd.color(f"0x{binder_hex}", f"{obj} and chain {bin_ch}")
    cmd.set("cartoon_flat_sheets", 0)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("cartoon_loop_radius", 0.4)

    # interface residues (< 5 A heavy-atom contact across chains)
    cmd.select("ifA", f"(({veg_sel}) within 5 of chain {bin_ch})")
    cmd.select("ifB", f"(chain {bin_ch} within 5 of ({veg_sel}))")
    cmd.show("sticks", "ifA"); cmd.show("spheres", "ifA")
    cmd.color("red", "ifA")
    cmd.set("sphere_scale", 0.45, "ifA"); cmd.set("stick_radius", 0.18, "ifA")
    cmd.show("sticks", "ifB"); cmd.show("spheres", "ifB")
    cmd.color("orange", "ifB")
    cmd.set("sphere_scale", 0.45, "ifB"); cmd.set("stick_radius", 0.18, "ifB")
    cmd.show("cartoon", "ifA"); cmd.set("cartoon_color", "grey80", "ifA")
    cmd.show("cartoon", "ifB"); cmd.set("cartoon_color", f"0x{binder_hex}", "ifB")

    # lighting / background
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)   # force opaque white background in ray+png
    cmd.set("depth_cue", 1); cmd.set("fog_start", 0.35)
    cmd.set("light_count", 3); cmd.set("specular", 0.45); cmd.set("shininess", 40)
    cmd.set("ambient", 0.30); cmd.set("direct", 0.70)
    cmd.set("ray_shadows", 1); cmd.set("antialias", 2); cmd.set("ray_trace_mode", 1)

    # orient & zoom on the binding interface
    cmd.select("iface_union", "ifA or ifB")
    cmd.orient("iface_union")
    cmd.zoom("iface_union", 18)
    cmd.turn("y", azim); cmd.turn("x", elev)

    cmd.ray(1400, 1400)
    cmd.png(out_png, 1400, 1400, dpi=300, ray=1)
    print("saved", out_png)
    cmd.reinitialize()


if __name__ == "__main__":
    print("== V1 (weak contact, monomer) ==")
    render_complex(V1_PDB, "2a6f97",
                   os.path.join(OUT, "GA_V1_contact.png"), azim=140, elev=18)
    print("== V4 (bridges dimer) ==")
    render_complex(V4_CIF, "2a9d8f",
                   os.path.join(OUT, "GA_V4_bridge.png"), azim=160, elev=20)
