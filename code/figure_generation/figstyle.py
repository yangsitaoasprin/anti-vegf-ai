# -*- coding: utf-8 -*-
"""Enforce a printable minimum text size on matplotlib figures.

WHY THIS EXISTS
---------------
Elsevier, "Artwork sizing":
    the lettering on the artwork should have a finished, printed size of 7 pt
    for normal text and no smaller than 6 pt for subscript and superscript
    characters. Smaller lettering will yield text that is hardly legible.

matplotlib sizes text in points relative to the FIGURE, not to the printed
page.  A figure authored W_auth inches wide and then placed at PRINT_W inches
prints its text at

    printed_pt = fontsize_pt * PRINT_W / W_auth

Every figure in this manuscript was authored far wider than it is placed, so
all of its text shrank -- by ~1.2x in the best case (Fig. 8) and ~2.4x in the
worst (Fig. S1).  The type was legible on screen and illegible on paper.

WHY IT CALIBRATES ITSELF
------------------------
W_auth is not figsize: rcParams sets savefig.bbox="tight", so the realised
raster is cropped to the drawn content and its width is only known after
rendering.  We therefore render, measure, rescale, and repeat until the
measured width stops moving.

WHY IT REPORTS OVERLAPS
-----------------------
Enlarging type is exactly the operation that creates collisions, so the helper
records text-vs-text intersections BEFORE and AFTER and reports the delta.
Only NEW overlaps indicate damage; pre-existing ones are the author's design.
"""
import io
import json
import os
import re
import tempfile

from matplotlib.text import Text

MIN_PT = 7.2          # target: 7 pt floor plus a little margin for rounding
MAX_ITER = 8
GROWTH_CAP = 1.35     # canvas growth above which the lift is declared divergent
_BUILD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "build_docx.py")

# Placed widths, mirroring build_docx.py's FIG_WIDTH.  Used ONLY when
# build_docx.py is not co-located -- e.g. inside the published `repro/`
# bundle, which ships the figure code without the manuscript builder.
# build_docx.py remains the single source of truth: whenever it is present
# it wins, so a width change there still propagates automatically.
_FALLBACK_WIDTH = {
    "Fig1_pipeline.png": 6.0,
    "Fig7_binding_pose.png": 6.5,
    "Fig8_interface_stats.png": 6.0,
    "Fig5_developability.png": 5.6,
    "Fig2_redesign_dG.png": 6.0,
    "Fig4_crossval.png": 6.0,
    "Fig3_full_length.png": 5.6,
    "Fig6_dimer.png": 5.6,
    "Graphical_Abstract.png": 6.2,
    "FigS1_apr_recovery.png": 6.5,
}


def _text_box(t, renderer):
    """Text-only bounding box, in display px.

    Annotation.get_window_extent() returns the UNION of the text box and the
    arrow patch.  An arrow's bbox is the large empty rectangle spanned by its
    two end points, so using it as "text area" reports collisions with labels
    the arrow merely flies over.  Call the base-class method to get the glyphs.
    """
    try:
        return Text.get_window_extent(t, renderer=renderer)
    except Exception:
        return t.get_window_extent(renderer=renderer)


def placed_widths():
    """{png basename: placed width in inches}, read from build_docx.py.

    build_docx.py is the single source of truth: if a figure's placed width
    changes there, this helper follows automatically and every figure is
    re-scaled against the new width on the next render.

    When build_docx.py is not co-located (the published `repro/` bundle ships
    the figure code without the manuscript builder) this falls back to
    _FALLBACK_WIDTH so the scripts still run standalone."""
    if not os.path.exists(_BUILD):
        return dict(_FALLBACK_WIDTH)
    txt = io.open(_BUILD, encoding="utf-8").read()
    files = eval(re.search(r"FIG_FILES\s*=\s*(\{[^}]*\})", txt).group(1))
    width = eval(re.search(r"FIG_WIDTH\s*=\s*(\{[^}]*\})", txt).group(1))
    return {files[k]: float(width[k]) for k in files}


def _texts(fig):
    out = []
    for t in fig.findobj(match=Text):
        try:
            s = t.get_text()
        except Exception:
            continue
        if not s or not s.strip() or not t.get_visible():
            continue
        try:
            if t.get_fontsize() is None:
                continue
        except Exception:
            continue
        out.append(t)
    return out


def _scale(fig, k):
    n = 0
    for t in _texts(fig):
        try:
            t.set_fontsize(float(t.get_fontsize()) * k)
            n += 1
        except Exception:
            pass
    return n


def _lift_to_floor(fig, floor):
    """Raise every text smaller than `floor` up to it; leave larger text alone.

    Lifting only small type (rather than scaling everything) is what makes this
    converge.  Scaling uniformly DIVERGES on figures whose long centred titles
    already reach the canvas edge: enlarging a title widens the tight-bbox
    canvas in proportion, so printed_pt = fontsize * PRINT_W / W_auth barely
    moves.  Measured on Fig. 4: a uniform x1.56 grew the canvas 1.29x and bought
    only 1.21x of printed size, then kept diverging.  Clamping leaves the big
    type -- the artists that set the width -- untouched."""
    n = 0
    for t in _texts(fig):
        try:
            s = float(t.get_fontsize())
            if s < floor - 1e-9:
                t.set_fontsize(floor)
                n += 1
        except Exception:
            pass
    return n


def _measure(fig, dpi):
    """Realised content width in inches (what bbox_inches='tight' produces).

    Returns None if the renderer refuses -- which is how divergence surfaces:
    once the loop has lifted type to a few hundred points, FreeType answers
    "raster overflow" instead of drawing a glyph.  Callers must handle None.
    """
    fd, p = tempfile.mkstemp(suffix=".png", prefix="_figstyle_")
    os.close(fd)
    try:
        fig.savefig(p, dpi=dpi, bbox_inches="tight")
        with open(p, "rb") as fh:
            hdr = fh.read(26)
        import struct
        w = struct.unpack(">II", hdr[16:24])[0]
    except Exception:
        return None
    finally:
        try:
            os.remove(p)
        except OSError:
            pass
    return w / float(dpi)


def overlaps(fig, min_frac=0.12):
    """Text-vs-text collisions: [(a, b, overlap_area_px, frac_of_smaller)]."""
    try:
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
    except Exception:
        return []
    items = []
    for t in _texts(fig):
        try:
            bb = _text_box(t, r)
        except Exception:
            continue
        if bb.width <= 0 or bb.height <= 0:
            continue
        items.append((t.get_text().strip()[:44].replace("\n", " / "), bb))
    out = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            si, bi = items[i]
            sj, bj = items[j]
            ow = min(bi.x1, bj.x1) - max(bi.x0, bj.x0)
            oh = min(bi.y1, bj.y1) - max(bi.y0, bj.y0)
            if ow <= 0 or oh <= 0:
                continue
            area = ow * oh
            small = min(bi.width * bi.height, bj.width * bj.height)
            frac = area / small if small else 0.0
            if frac >= min_frac:
                out.append((si, sj, int(area), round(frac, 3)))
    return out


def normalize(fig, png_name, out_dir=None, dpi=300, min_pt=MIN_PT):
    """Scale every Text in `fig` so nothing prints below `min_pt` pt.

    Call immediately before saving.  Returns a report dict; when `out_dir` is
    given, also writes <png stem>_typography.json into it so the submission gate
    can assert the floor after the fact."""
    w_in = placed_widths().get(png_name)
    rep = {"png": png_name, "placed_w_in": w_in, "target_min_pt": min_pt,
           "status": "ok", "iterations": [], "new_overlaps": []}
    if not w_in:
        rep["status"] = "no placed width for %s" % png_name
        return rep

    before = overlaps(fig)
    wauth0 = None
    for it in range(MAX_ITER):
        wauth = _measure(fig, dpi)
        if wauth is None:
            rep["status"] = ("diverged at iteration %d (renderer refused; the "
                             "lift outgrew the canvas)" % it)
            break
        if wauth0 is None:
            wauth0 = wauth
        elif wauth > wauth0 * GROWTH_CAP:
            # Lifting raised the canvas width, which raised the floor, which
            # lifts again.  This figure cannot be fixed by re-scaling type: it
            # carries too much text for the width it is placed at.
            rep["status"] = ("diverged: canvas grew %.2fx to %.2f in"
                             % (wauth / wauth0, wauth))
            rep["divergence"] = {"wauth0_in": round(wauth0, 3),
                                 "wauth_in": round(wauth, 3),
                                 "growth": round(wauth / wauth0, 3),
                                 "placed_w_in": w_in,
                                 "min_fontsize_pt": round(
                                     min(float(t.get_fontsize()) for t in _texts(fig)), 2),
                                 "required_floor_pt": round(
                                     min_pt * wauth / w_in, 2)}
            break
        sizes = [float(t.get_fontsize()) for t in _texts(fig)]
        if not sizes:
            rep["status"] = "no text"
            return rep
        smin = min(sizes)
        floor = min_pt * wauth / w_in          # pt that prints at exactly min_pt
        printed = smin * w_in / wauth
        rep["iterations"].append({"wauth_in": round(wauth, 3),
                                  "canvas_pt_floor": round(floor, 2),
                                  "min_fontsize_pt": round(smin, 2),
                                  "printed_min_pt": round(printed, 3),
                                  "lifted": _lift_to_floor(fig, floor)})
        if printed >= min_pt - 0.005:
            break

    wauth_f = _measure(fig, dpi)
    sizes = [float(t.get_fontsize()) for t in _texts(fig)]
    if wauth_f is not None and sizes:
        rep["authoring_w_in"] = round(wauth_f, 3)
        rep["min_fontsize_pt"] = round(min(sizes), 3)
        rep["printed_min_pt"] = round(min(sizes) * w_in / wauth_f, 3)
        rep["max_fontsize_pt"] = round(max(sizes), 3)
    after = overlaps(fig)
    key = lambda o: tuple(sorted((o[0], o[1])))
    seen = {key(o) for o in before}
    rep["new_overlaps"] = [o for o in after if key(o) not in seen]
    rep["pre_existing_overlaps"] = len(before)
    # The full post-lift set, not just the delta. "New" answers "did enlarging
    # the type do damage?"; for a submission we additionally require that no
    # label touches another at all, so the gate reads this field.
    rep["overlaps_after"] = after
    if rep["status"] == "ok":
        if rep.get("printed_min_pt", 0) < min_pt - 0.005:
            rep["status"] = "did not converge"
        if rep["new_overlaps"]:
            rep["status"] = "new text overlaps introduced"

    if out_dir:
        path = os.path.join(out_dir,
                            os.path.splitext(png_name)[0] + "_typography.json")
        with io.open(path, "w", encoding="utf-8") as fh:
            fh.write(json.dumps(rep, indent=2, ensure_ascii=False))
        rep["sidecar"] = path
    return rep


def report_line(rep):
    return ("%-24s placed %.2f in | authoring %.2f in | %d->%.1f pt | "
            "printed min %.2f pt | %s"
            % (rep.get("png", "?"), rep.get("placed_w_in") or 0,
               rep.get("authoring_w_in") or 0, rep.get("min_fontsize_pt") or 0,
               rep.get("max_fontsize_pt") or 0, rep.get("printed_min_pt") or 0,
               rep.get("status", "?")))
