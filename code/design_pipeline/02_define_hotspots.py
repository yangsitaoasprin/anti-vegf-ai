# =============================================================================
# Author : Yang Sitao (杨四涛)
# Affiliation : School of Pharmaceutical Sciences, Dali University;
#               The Third Affiliated Hospital of Dali University
# Contact : yangsitaoasprin@126.com
# License : MIT (see ../../LICENSE)
# =============================================================================
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
02_define_hotspots.py
从 VEGF-A / 受体复合物结构中，通过界面残基分析定义 RFdiffusion 的 hotspot 残基。

原理:
    界面残基 = 配体链(VEGF-A)中，任一重原子与受体链任一重原子距离 < cutoff 的残基。
    这些残基是 VEGF-A 与受体结合的"热点"，作为设计新结合物的靶向位点。

用法:
    python 02_define_hotspots.py --complex ../data/1FLT.pdb \
        --ligand-chain W --receptor-chain X --cutoff 5.0 --out ../data/hotspots.txt

依赖:
    Biopython (pip install biopython)。未安装时回退到内置纯文本解析器。

严谨性:
    - 输出 hotspot 及其到受体的最小距离，可复核。
    - 同时输出一份"打乱/随机 hotspot"作为阴性对照（供 H1 假设检验）。
"""
import argparse
import os
import random
import sys
from math import sqrt


def parse_atoms_plain(path, chains):
    """纯文本解析 PDB，返回 {chain: [(resid, resname, atomname, x, y, z), ...]}。"""
    data = {c: [] for c in chains}
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            ch = line[21]
            if ch not in chains:
                continue
            try:
                resid = int(line[22:26])
                resname = line[17:20].strip()
                atom = line[12:16].strip()
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            if atom.startswith("H"):  # 忽略氢
                continue
            data[ch].append((resid, resname, atom, x, y, z))
    return data


def interface_residues(lig_atoms, rec_atoms, cutoff):
    """返回配体链界面残基: {(resid, resname): min_distance}。O(N*M) 朴素实现，靶点规模足够。"""
    cutoff2 = cutoff * cutoff
    rec_coords = [(a[3], a[4], a[5]) for a in rec_atoms]
    hits = {}
    for resid, resname, atom, x, y, z in lig_atoms:
        best = None
        for rx, ry, rz in rec_coords:
            d2 = (x - rx) ** 2 + (y - ry) ** 2 + (z - rz) ** 2
            if d2 < cutoff2:
                d = sqrt(d2)
                if best is None or d < best:
                    best = d
        if best is not None:
            key = (resid, resname)
            if key not in hits or best < hits[key]:
                hits[key] = best
    return hits


def main():
    ap = argparse.ArgumentParser(description="定义 VEGF-A 界面 hotspot 残基")
    ap.add_argument("--complex", required=True, help="复合物 PDB 文件")
    ap.add_argument("--ligand-chain", required=True, help="VEGF-A 链 (配体)")
    ap.add_argument("--receptor-chain", required=True, help="受体链")
    ap.add_argument("--cutoff", type=float, default=5.0, help="界面距离阈值 (Å)")
    ap.add_argument("--out", default="../data/hotspots.txt")
    args = ap.parse_args()

    if not os.path.exists(args.complex):
        sys.exit(f"[错误] 找不到结构文件: {args.complex}")

    chains = [args.ligand_chain, args.receptor_chain]
    # 优先 Biopython；失败回退纯文本
    try:
        from Bio.PDB import PDBParser, NeighborSearch  # noqa
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("cplx", args.complex)
        model = structure[0]
        lig = [a for a in model[args.ligand_chain].get_atoms() if a.element != "H"]
        rec = [a for a in model[args.receptor_chain].get_atoms() if a.element != "H"]
        ns = NeighborSearch(rec)
        hits = {}
        for atom in lig:
            near = ns.search(atom.coord, args.cutoff)
            if near:
                res = atom.get_parent()
                resid = res.id[1]
                resname = res.resname
                d = min(float(((atom.coord - n.coord) ** 2).sum()) ** 0.5 for n in near)
                key = (resid, resname)
                if key not in hits or d < hits[key]:
                    hits[key] = d
        print("[信息] 使用 Biopython 完成界面分析。")
    except ImportError:
        print("[信息] 未安装 Biopython，使用内置纯文本解析器。")
        data = parse_atoms_plain(args.complex, chains)
        hits = interface_residues(data[args.ligand_chain], data[args.receptor_chain], args.cutoff)

    if not hits:
        sys.exit("[错误] 未发现界面残基，请核对链标识与 cutoff。")

    hotspots = sorted(hits.keys(), key=lambda k: k[0])
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(f"# VEGF-A hotspot 残基 (链 {args.ligand_chain}, cutoff {args.cutoff}Å)\n")
        fh.write(f"# 格式: 链{args.ligand_chain}残基号  残基名  到受体最小距离(Å)\n")
        for (resid, resname) in hotspots:
            fh.write(f"{args.ligand_chain}{resid}\t{resname}\t{hits[(resid, resname)]:.2f}\n")
    print(f"[成功] 发现 {len(hotspots)} 个界面残基 -> {args.out}")
    # RFdiffusion 风格 hotspot 串
    hs_str = ",".join(f"{args.ligand_chain}{r}" for r, _ in hotspots)
    print(f"[RFdiffusion hotspot] {hs_str}")

    # 阴性对照：随机残基（同数量），供 H1 假设检验
    all_res = sorted({r for r, _ in [(a[0], a[1]) for a in []]})  # placeholder
    ctrl_path = args.out.replace(".txt", "_negctrl.txt")
    random.seed(42)
    # 从配体链全部残基中随机抽同等数量（排除真实 hotspot）
    try:
        data = parse_atoms_plain(args.complex, [args.ligand_chain])
        all_ids = sorted({a[0] for a in data[args.ligand_chain]})
        hot_ids = {r for r, _ in hotspots}
        pool = [i for i in all_ids if i not in hot_ids]
        ctrl = random.sample(pool, min(len(hotspots), len(pool)))
        with open(ctrl_path, "w", encoding="utf-8") as fh:
            fh.write("# 阴性对照 hotspot（随机残基，用于 H1 假设检验）\n")
            for i in sorted(ctrl):
                fh.write(f"{args.ligand_chain}{i}\n")
        print(f"[对照] 随机阴性对照 hotspot -> {ctrl_path}")
    except Exception as e:
        print(f"[警告] 生成阴性对照失败: {e}")


if __name__ == "__main__":
    main()
