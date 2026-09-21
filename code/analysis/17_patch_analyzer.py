#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
17_patch_analyzer.py
====================
用 Schrödinger BioLuminate Patch Analyzer 给候选微蛋白补「结构维度」的
聚集风险评估 —— 这是 WALTZ / APR-Score 两条序列线都给不了的信息：

    高风险六肽到底露在表面，还是埋在疏水核里？

为什么需要它
------------
WALTZ 与 APR-Score 都是**纯序列**预测器：它们说「这段序列有聚集倾向」，
但不问这段序列在折叠后的蛋白里能不能被看见。一个埋在核心里的 APR 和
一个摆在表面的 APR，可开发性后果完全不同。Patch Analyzer 补的正是这一层。

它同时带来两个**血统不同**的序列聚集预测器（与 WALTZ-DB、APR-Score 均不同源）：
    * Zyggregator  (Tartaglia / Chiti 系)
    * Aggrescan    (Conchillo-Sole 系)
以及结构量：每残基 SASA、侧链 SASA、accessibility(相对暴露%)、
表面 patch（疏水/正电/负电，面积 Å²）与 AggScore。

两个状态
--------
    mono : binder 单链（从复合物中取出，坐标不变）
           -> 采信 Zyggregator / Aggrescan / SASA / patches
           -> 这是「游离态」，可开发性相关
    cpx  : binder + 靶点复合物
           -> 只采信 SASA / accessibility，用于回答「靶点遮住了多少」
           -> 不采信 cpx 里的 Zyggregator：该算法带净电荷项，
              Schr 实现用的是整条 CT 的总电荷（两条链之和），非 binder 自有电荷

用法
----
    python 17_patch_analyzer.py                # 全部跑
    python 17_patch_analyzer.py --only V1 V2   # 只跑部分
    python 17_patch_analyzer.py --skip-prep    # 跳过 prepwizard（已有 maegz）

依赖：Schrödinger 2025-2（prepwizard + BioLuminate Patch Analyzer）
"""
import argparse
import json
import os
import subprocess
import sys
import io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)

SCHRODINGER = r"C:\Program Files\Schrodinger2025-2"
PREPWIZARD = os.path.join(SCHRODINGER, "utilities", "prepwizard.exe")
SCHPY = os.path.join(SCHRODINGER, "internal", "bin", "python.exe")

SRC_DIR = os.path.join(PROJ, "designs", "candidates_6", "wetlab", "structures")
WORK = os.path.join(PROJ, "results", "patch_analyzer")
OUT_JSON = os.path.join(WORK, "patch_analyzer_raw.json")

# 靶点标签：VEGF-A RBD 的 N 端特征片段（设计中使用的 TARGET 前段）
VEGF_TAG = "HEVVKFMDVYQRSYCH"

# ---------------------------------------------------------------
# 候选 -> 复合物结构
# ---------------------------------------------------------------
CONSTRUCTS = [
    # tag,  面板名,                  复合物文件,                              角色
    ("V1", "#1_T0.3_s132",           "1star_fullVEGF_1T0.3_s132_complex.pdb", "1* 去A重设计"),
    ("V2", "#2c_T0.1_s18",           "2star_fullVEGF_2c_T0.1_s18_complex.pdb", "2* 保守微调"),
    ("V3", "vegf_len60_1_T0.2_s1",   "vegf_len60_1_T0.2_s1.complex.pdb",      "3  未修饰"),
    ("V4", "T0.1_s35",               "4star_fullVEGF_T0.1_s35_complex.pdb",   "4* 去淀粉样"),
    ("V5", "vegf_len100_3_T0.1_s2",  "vegf_len100_3_T0.1_s2.complex.pdb",     "5  未修饰"),
    ("P1", "vegf_len80_1_T0.2_s2",   "vegf_len80_1_T0.2_s2.complex.pdb",      "父本 #1"),
    ("P2", "vegf_len100_4_T0.2_s1",  "vegf_len100_4_T0.2_s1.complex.pdb",     "父本 #2"),
    ("P4", "vegf_len80_3_T0.1_s2",   "vegf_len80_3_T0.1_s2.complex.pdb",      "父本 #4"),
]

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "HID": "H", "HIE": "H", "HIP": "H",
}


_HOME = os.environ.get("USERPROFILE") or os.path.expanduser("~")
_HOMEPATH = os.environ.get("HOMEPATH") or os.path.splitdrive(_HOME)[1] or "\\Users\\Public"


def clean_env():
    """最小化环境：绕开宿主会话注入的 sitecustomize.py 与超长环境变量。"""
    return {
        "PATH": (r"C:\Windows\System32;C:\Windows;" + SCHRODINGER + ";"
                 + os.path.join(SCHRODINGER, "utilities") + ";"
                 + os.path.join(SCHRODINGER, "mmshare-v7.0", "bin", "Windows-x64")),
        "SCHRODINGER": SCHRODINGER,
        "SCHRODINGER_PYTHON": SCHPY,
        "LM_LICENSE_FILE": os.path.join(SCHRODINGER, "licenses", "25-2.lic"),
        "PYTHONNOUSERSITE": "1",
        "SYSTEMROOT": r"C:\Windows",
        "SystemRoot": r"C:\Windows",
        "COMSPEC": r"C:\Windows\System32\cmd.exe",
        "TEMP": os.environ.get("TEMP", r"C:\Windows\Temp"),
        "TMP": os.environ.get("TMP", r"C:\Windows\Temp"),
        "NUMBER_OF_PROCESSORS": os.environ.get("NUMBER_OF_PROCESSORS", "8"),
        "PROCESSOR_ARCHITECTURE": "AMD64",
        "OS": "Windows_NT",
        "USERPROFILE": _HOME,
        "HOMEDRIVE": "C:",
        "HOMEPATH": _HOMEPATH,
        "USERNAME": os.environ.get("USERNAME", "user"),
        "WINDIR": r"C:\Windows",
    }


# ---------------------------------------------------------------
# 1. 拆链 / 找 binder
# ---------------------------------------------------------------
def read_chains(pdb_path):
    """
    返回 {chain: {"seq": 一字母序列, "lines": 全部 ATOM 行, "nres": 残基数}}

    ⚠ 去重只用于拼序列；写出结构时必须是**全部** ATOM 行。
      曾在此踩坑：把去重后的行当结构写出 -> 每残基只剩 1 个原子，
      prepwizard 产出 4 原子/残基的残废结构，Patch Analyzer 仍能跑但
      所有数值无意义（且 get_connected_components 退化成 80 个单残基组分，
      Zyggregator 因 len<5 全部返回 None）。
    """
    chains = {}
    with open(pdb_path, "r", errors="replace") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[16] not in (" ", "A"):      # altLoc
                continue
            ch = line[21]
            d = chains.setdefault(ch, {"seq": [], "lines": [], "seen": set()})
            d["lines"].append(line)
            rkey = line[22:27]                   # resSeq + iCode
            if rkey not in d["seen"]:
                d["seen"].add(rkey)
                d["seq"].append(AA3TO1.get(line[17:20].strip(), "X"))
    for c, d in chains.items():
        d["seq"] = "".join(d["seq"])
        d["nres"] = len(d["seen"])
        d["natm"] = len(d["lines"])
    return chains


def pick_binder(chains):
    """binder = 不含 VEGF 特征片段的那条链。"""
    cand = []
    for ch, d in chains.items():
        if VEGF_TAG in d["seq"]:
            continue
        cand.append((d["nres"], ch, d["seq"]))
    if len(cand) != 1:
        raise RuntimeError("binder 链判定失败: %r"
                           % [(c, chains[c]["nres"]) for _, c, _ in cand])
    return cand[0][1], cand[0][2]


# ---------------------------------------------------------------
# 2. 内层 Schrödinger 脚本
# ---------------------------------------------------------------
INNER = r'''
# -*- coding: utf-8 -*-
import json, sys
from schrodinger import structure
from schrodinger.application.bioluminate.patch_utils import patch_finder as pf
from schrodinger.application.bioluminate.patch_utils import settings as pset

basename, outjson = sys.argv[1], sys.argv[2]

pre = pf.PreAnalyzer.readAndRun(basename)
analysis = pre.getAnalysis()
finder = pf.PatchFinder(analysis)
Patches, res_data = finder.calculate(pset.PatchSettings())

# --- 完整残基集合：直接向 _calcAggregation 要，覆盖"表面无顶点"的完全埋藏残基 ---
# （finder._res_data_by_vertex 只含至少占到一个表面顶点的残基；完全埋藏者会缺失，
#   而"APR 是否被埋在核心"正是本项目要回答的问题，不能漏。）
all_seq, all_struc = pf.get_connected_components(pre._struc)
agg_by_res = pre._calcAggregation(all_seq, all_struc)

all_res = {}
for res, info in agg_by_res.items():
    all_res[(info.chain, info.resnum, info.inscode)] = {
        "chain": info.chain, "resnum": info.resnum,
        "inscode": info.inscode, "resname": info.resname.strip(),
        "zyggregator": info.zyggregator, "aggrescan": info.aggrescan,
        "sasa": 0.0, "side_chain_sasa": 0.0, "accessibility": 0.0,
        "on_surface": False, "aggscore": 0.0,
        "topography": info.topography,
        "reactive_type": (info.reactive_type or None),
    }

n_on_surface = 0
for r in finder._res_data_by_vertex:
    if r is None:
        continue
    key = (r.chain, r.resnum, r.inscode)
    if key in all_res:
        d = all_res[key]
        if not d["on_surface"]:
            n_on_surface += 1
        d["sasa"] = r.sasa
        d["side_chain_sasa"] = r.side_chain_sasa
        d["accessibility"] = (r.accessibility if r.accessibility is not None else 0.0)
        d["on_surface"] = True

# --- patch 残基：补 aggscore ---
for rd in res_data:
    key = (rd.chain, rd.resnum, rd.inscode)
    if key in all_res:
        all_res[key]["aggscore"] = rd.aggscore
    else:
        all_res[key] = {
            "chain": rd.chain, "resnum": rd.resnum,
            "inscode": rd.inscode, "resname": rd.resname.strip(),
            "zyggregator": rd.zyggregator, "aggrescan": rd.aggrescan,
            "sasa": rd.sasa, "side_chain_sasa": rd.side_chain_sasa,
            "accessibility": (rd.accessibility if rd.accessibility is not None else 0.0),
            "on_surface": True, "aggscore": rd.aggscore,
            "topography": rd.topography,
            "reactive_type": (rd.reactive_type or None),
        }

patches_out = []
for p in Patches:
    try:
        resnames = [r.resname.strip() for r in p.residues]
        contrib = [float(c) for c in p.contribution]
    except Exception:
        resnames, contrib = [], []
    patches_out.append({
        "type": p.type.long_name,
        "size": float(p.size),
        "n_vertices": len(p.vertices),
        "smoothed_mean": (sum(p.vtxvals) / len(p.vtxvals)) if p.vtxvals else None,
        "smoothed_max": (max(p.vtxvals) if p.vtxvals else None),
        "energy": (float(p.energy) if p.energy is not None else None),
        "residues": resnames,
        "contribution": contrib,
    })

props = {}
for k, v in vars(finder.prot_properties).items():
    try:
        json.dumps(v)
        props[k] = v
    except Exception:
        props[k] = repr(v)

out = {
    "basename": basename,
    "n_components": len(all_seq),
    "component_lengths": [len(s) for s in all_seq],
    "n_res_total": len(all_res),
    "n_on_surface": n_on_surface,
    "residues": sorted(all_res.values(), key=lambda d: (d["chain"], d["resnum"])),
    "patches": patches_out,
    "properties": props,
}
with open(outjson, "w", encoding="utf-8") as fh:
    json.dump(out, fh, ensure_ascii=False)
sys.stdout.write("OK %s residues=%d/%d-onSurf components=%d patches=%d\n"
                 % (basename, len(all_res), n_on_surface, len(all_seq), len(Patches)))
'''


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--skip-prep", action="store_true")
    a = ap.parse_args()

    os.makedirs(WORK, exist_ok=True)
    inner_py = os.path.join(WORK, "_inner_patch.py")
    with open(inner_py, "w", encoding="utf-8") as fh:
        fh.write(INNER)

    todo = [c for c in CONSTRUCTS if (a.only is None or c[0] in a.only)]
    env = clean_env()
    result = {}

    for tag, name, fname, role in todo:
        src = os.path.join(SRC_DIR, fname)
        if not os.path.exists(src):
            print("[跳过] %s 结构缺失: %s" % (tag, fname), flush=True)
            continue
        print("=" * 78, flush=True)
        print("%s  %s  (%s)" % (tag, name, role), flush=True)

        chains = read_chains(src)
        bchain, bseq = pick_binder(chains)
        tchains = [c for c in chains if c != bchain]
        print("  链: %s   binder=%s(%d aa)   target=%s"
              % ({c: "%dres/%datm" % (chains[c]["nres"], chains[c]["natm"])
                  for c in chains},
                 bchain, len(bseq), tchains), flush=True)
        # 完整性自检：重原子 PDB 应 >= 4 原子/残基（否则说明结构被写残了）
        for c, d in chains.items():
            apr = d["natm"] / max(d["nres"], 1)
            if apr < 4.0:
                raise RuntimeError(
                    "%s 链 %s 仅 %.2f 原子/残基，结构不完整，终止" % (tag, c, apr))

        # --- 写两种输入 ---
        base_mono = os.path.join(WORK, "%s_mono" % tag)
        base_cpx = os.path.join(WORK, "%s_cpx" % tag)
        with open(base_mono + ".pdb", "w", encoding="utf-8") as fh:
            fh.writelines(chains[bchain]["lines"])
            fh.write("END\n")
        with open(base_cpx + ".pdb", "w", encoding="utf-8") as fh:
            for c in sorted(chains):
                fh.writelines(chains[c]["lines"])
            fh.write("END\n")
        # ASL —— 必须与传给 PreAnalyzer.readAndRun 的 basename 同名
        # （basename = <base>_prepped，故写 <base>_prepped.txt）
        with open(base_mono + "_prepped.txt", "w", encoding="utf-8") as fh:
            fh.write('chain.name "%s"' % bchain)
        with open(base_cpx + "_prepped.txt", "w", encoding="utf-8") as fh:
            fh.write(" OR ".join('chain.name "%s"' % c for c in sorted(chains)))

        rec = {"tag": tag, "name": name, "role": role,
               "binder_chain": bchain, "binder_seq": bseq,
               "target_chains": tchains,
               "target_len": {c: chains[c]["nres"] for c in tchains},
               "structure": fname}

        for state, base in (("mono", base_mono), ("cpx", base_cpx)):
            mae = base + "_prepped.maegz"
            if (not a.skip_prep) or (not os.path.exists(mae)):
                r = subprocess.run(
                    [PREPWIZARD, base + ".pdb", mae,
                     "-fillsidechains", "-disulfides", "-propka_pH", "7.4",
                     "-HOST", "localhost", "-WAIT"],
                    env=env, capture_output=True, text=True, cwd=WORK, timeout=3600)
                if r.returncode != 0 or not os.path.exists(mae):
                    print("  [失败] prepwizard %s rc=%d" % (state, r.returncode), flush=True)
                    print("         " + (r.stderr or "")[-600:].replace("\n", "\n         "), flush=True)
                    continue
                print("  prepwizard %s  OK" % state, flush=True)
            else:
                print("  prepwizard %s  复用" % state, flush=True)

            # maegz + ASL -> 交给内层 Schr python
            js = base + "_patch.json"
            r = subprocess.run(
                [SCHPY, "-E", "-s", inner_py, base + "_prepped", js],
                env=env, capture_output=True, text=True, cwd=WORK, timeout=3600)
            if r.returncode != 0 or not os.path.exists(js):
                print("  [失败] PatchAnalyzer %s rc=%d" % (state, r.returncode), flush=True)
                print("         " + (r.stdout or "")[-400:].replace("\n", "\n         "), flush=True)
                print("         " + (r.stderr or "")[-900:].replace("\n", "\n         "), flush=True)
                continue
            print("  " + (r.stdout or "").strip().splitlines()[-1], flush=True)
            with open(js, "r", encoding="utf-8") as fh:
                rec[state] = json.load(fh)

        result[tag] = rec

    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump(result, fh, ensure_ascii=False)
    print("\n写入 " + OUT_JSON, flush=True)


if __name__ == "__main__":
    main()
