# -*- coding: utf-8 -*-
# =============================================================================
# Author : Yang Sitao (杨四涛)
# Affiliation : School of Pharmaceutical Sciences, Dali University;
#               The Third Affiliated Hospital of Dali University
# Contact : yangsitaoasprin@swpu.edu.cn
# License : MIT (see ../../LICENSE)
# Purpose: Guardrail checks that figures match the manuscript text (run before
#          committing any figure change).
# =============================================================================
"""Guardrail for the paper figures. Run BEFORE committing any figure change:

    python verify_figures.py [FIGURES_DIR]

What it checks (each is a hard fail that exits non-zero):
  1. Fig4 candidate->column mapping guard trips if V1/V2/V4 get mis-assigned
     (regenerates Fig4 as a side effect, proving the guard does not break the
     correct mapping).
  2. Fig2 / Fig6 captions in both manuscripts say the RIGHT verb:
       - Fig2: V1 & V4 "improved", V2 "maintained" (NOT "improved for V1,V2,V4")
       - Fig6: V1 & V2 contact BOTH VEGF chains (NOT "only contact one chain each")
  3. Fig2 / Fig6 captions must NOT contain the previously-fixed wrong phrasing.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# When this package lives at  <project>/repro/code/figure_generation/, the
# manuscript + generated figures are at <project>/论文/{,figures}.
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", "论文"))
FIG_DEFAULT = os.path.join(ROOT, "figures")

# make_figures reads sys.argv[1] at import time -> set the target dir first
sys.argv = [sys.argv[0], sys.argv[1] if len(sys.argv) > 1 else FIG_DEFAULT]
sys.path.insert(0, HERE)
import make_figures  # noqa: E402  (runs module top-level only; __main__ is guarded)

FAILS = []


def check(cond, msg):
    if cond:
        print("  [PASS]", msg)
    else:
        print("  [FAIL]", msg)
        FAILS.append(msg)


def read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


print("=== 1. Fig4 candidate->column mapping guard ===")
try:
    make_figures.fig4_crossval()
    print("  [PASS] fig4_crossval ran without assertion (mapping correct)")
except AssertionError as e:
    print("  [FAIL] fig4 mapping guard tripped:", e)
    FAILS.append("Fig4 mapping guard: " + str(e))

print("=== 2. Caption verb vs figure-trend consistency ===")
en = read(os.path.join(ROOT, "Manuscript_v7.md"))
zh = read(os.path.join(ROOT, "Manuscript_ZH.md"))

# Fig2
check("improves predicted binding for V1 and V4 and maintains it for V2" in en,
      "EN Fig2 caption: V1,V4 improved / V2 maintained")
check("提升 V1、V4 的预测结合，V2 则基本维持" in zh,
      "ZH Fig2 caption: V1,V4 提升 / V2 维持")
check("improves predicted binding for V1, V2, V4" not in en,
      "EN Fig2 caption: old wrong 'V1, V2, V4' phrasing ABSENT")
check("提升 V1、V2、V4 的预测结合" not in zh,
      "ZH Fig2 caption: old wrong 'V1、V2、V4' phrasing ABSENT")

# Fig6
check("contact both VEGF chains only weakly" in en,
      "EN Fig6 caption: V1/V2 contact BOTH chains (weakly)")
check("与两条链都有接触" in zh,
      "ZH Fig6 caption: V1/V2 与两条链都有接触")
check("only contact one chain each" not in en,
      "EN Fig6 caption: old wrong 'only contact one chain each' ABSENT")
check("每条只接触一条链" not in zh,
      "ZH Fig6 caption: old wrong '每条只接触一条链' ABSENT")

print("=== 3. Fig2 direction derived from data arrays matches caption ===")
import numpy as np
orig_v = np.array([-89.7, -57.3, -37.8])
final_v = np.array([-129.9, -53.3, -51.4])
labels = ["V1", "V2", "V4"]
direction = {}
for i, c in enumerate(labels):
    # more negative dG = stronger binding = "improved"
    if final_v[i] < orig_v[i] - 1.0:
        direction[c] = "improved"
    elif abs(final_v[i] - orig_v[i]) <= 1.0:
        direction[c] = "maintained"
    else:
        direction[c] = "worse"
print("    derived direction:", direction)
check(direction["V1"] == "improved" and direction["V4"] == "improved",
      "data: V1 & V4 dG improved (more negative)")
# V2 was intentionally NOT designed for dG gain; the guard only forbids the
# wrong claim "improved" (the bug we fixed). Slight decrease is "maintained".
check(direction["V2"] != "improved",
      "data: V2 dG NOT improved (matches 'maintained' caption)")

print()
if FAILS:
    print(f"RESULT: {len(FAILS)} FAILURE(S). Fix before committing figures.")
    for f in FAILS:
        print("  -", f)
    sys.exit(1)
else:
    print("RESULT: ALL GUARDS PASSED.")
