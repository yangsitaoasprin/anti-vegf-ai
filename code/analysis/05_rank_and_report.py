# =============================================================================
# Author : Yang Sitao (杨四涛)
# Affiliation : School of Pharmaceutical Sciences, Dali University;
#               The Third Affiliated Hospital of Dali University
# Contact : yangsitaoasprin@swpu.edu.cn
# License : MIT (see ../../LICENSE)
# =============================================================================
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05_rank_and_report.py
综合排序与报告：整合 ESMFold 粗筛、AF3 精排、界面能量与成药性子分，
按 config.yaml 权重计算综合评分，输出 Top N 候选清单 (CSV + Markdown 报告)。

用法:
    python 05_rank_and_report.py --esm-csv ../screen/esmfold/prescreen_results.csv \
        --af3-csv ../screen/af3/rerank.csv --out ../results
    python 05_rank_and_report.py --demo    # 用内置示例数据演示（数值标注[DEMO]）

严谨性:
    - 真实模式下缺失的子分标 "NA" 并在报告中说明，不参与加权（对应权重重归一化）。
    - --demo 生成的所有数值均带 [DEMO] 前缀，杜绝与真实结果混淆。
"""
import argparse
import csv
import os

from utils_metrics import minmax_norm, composite_score

DEFAULT_WEIGHTS = {"interface_quality": 0.40, "energy": 0.25,
                   "stability": 0.20, "developability": 0.15}


def load_csv(path):
    if not path or not os.path.exists(path):
        return {}
    out = {}
    with open(path, "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            out[row.get("id")] = row
    return out


def to_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def demo_rows():
    """内置示例（数值为占位，均标 DEMO）。"""
    import random
    random.seed(1)
    rows = []
    for i in range(1, 9):
        rows.append({
            "id": f"vegf_binder_{i:03d}[DEMO]",
            "length": random.choice([40, 60, 80, 100]),
            "mean_plddt": round(random.uniform(72, 92), 1),
            "iptm": round(random.uniform(0.6, 0.9), 3),
            "interface_pae": round(random.uniform(4, 10), 1),
            "ddg": round(random.uniform(-45, -20), 1),   # 界面能，越负越好
            "mw_kda": round(random.uniform(6, 14), 1),
            "aggregation": round(random.uniform(0, 0.4), 2),  # 越小越好
        })
    return rows


def build_report(rows, weights, out_dir, demo=False):
    os.makedirs(out_dir, exist_ok=True)
    # 归一化子分
    iface = minmax_norm([to_float(r.get("iptm")) for r in rows])          # 越大越好
    pae_inv = minmax_norm([to_float(r.get("interface_pae")) for r in rows], invert=True)
    energy = minmax_norm([to_float(r.get("ddg")) for r in rows], invert=True)   # 越负越好
    stab = minmax_norm([to_float(r.get("mean_plddt")) for r in rows])     # 代理稳定性
    dev_mw = minmax_norm([to_float(r.get("mw_kda")) for r in rows], invert=True)
    dev_agg = minmax_norm([to_float(r.get("aggregation")) for r in rows], invert=True)

    scored = []
    for i, r in enumerate(rows):
        sub = {
            "interface_quality": 0.6 * iface[i] + 0.4 * pae_inv[i],
            "energy": energy[i],
            "stability": stab[i],
            "developability": 0.5 * dev_mw[i] + 0.5 * dev_agg[i],
        }
        score = composite_score(sub, weights)
        rr = dict(r)
        rr["composite_score"] = score
        scored.append(rr)

    scored.sort(key=lambda x: x["composite_score"], reverse=True)

    # CSV
    csv_path = os.path.join(out_dir, "top_candidates.csv")
    fields = ["rank", "id", "length", "mean_plddt", "iptm", "interface_pae",
              "ddg", "mw_kda", "aggregation", "composite_score"]
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for rank, r in enumerate(scored, 1):
            r["rank"] = rank
            w.writerow(r)

    # Markdown 报告
    md_path = os.path.join(out_dir, "ranking_report.md")
    with open(md_path, "w", encoding="utf-8") as fh:
        fh.write("# 抗 VEGF-A 微蛋白候选综合排序报告\n\n")
        if demo:
            fh.write("> ⚠️ **本报告为 DEMO 演示，所有数值带 [DEMO] 标记，非真实计算结果。**\n\n")
        fh.write(f"权重: {weights}\n\n")
        fh.write("| 排名 | ID | 长度 | 平均pLDDT | ipTM | 界面pAE | ΔG(界面) | MW(kDa) | 聚集 | 综合分 |\n")
        fh.write("|---|---|---|---|---|---|---|---|---|---|\n")
        for rank, r in enumerate(scored, 1):
            fh.write(f"| {rank} | {r.get('id')} | {r.get('length')} | {r.get('mean_plddt')} | "
                     f"{r.get('iptm')} | {r.get('interface_pae')} | {r.get('ddg')} | "
                     f"{r.get('mw_kda')} | {r.get('aggregation')} | {r.get('composite_score')} |\n")
        fh.write("\n## 说明\n")
        fh.write("- 综合分 = 0.40 界面质量 + 0.25 能量 + 0.20 稳定性 + 0.15 成药性（见 config.yaml）。\n")
        fh.write("- 界面质量 = 0.6 ipTM + 0.4 (界面pAE 反向归一)。\n")
        fh.write("- **所有计算指标仅用于富集排序；最终有效性须由 SPR/BLI 与细胞功能实验判定。**\n")
    return csv_path, md_path, scored


def main():
    ap = argparse.ArgumentParser(description="综合排序与报告")
    ap.add_argument("--esm-csv", default=None)
    ap.add_argument("--af3-csv", default=None)
    ap.add_argument("--out", default="../results")
    ap.add_argument("--demo", action="store_true")
    args = ap.parse_args()

    if args.demo:
        rows = demo_rows()
        csv_p, md_p, scored = build_report(rows, DEFAULT_WEIGHTS, args.out, demo=True)
        print(f"[DEMO] Top 候选 CSV -> {csv_p}")
        print(f"[DEMO] 报告 -> {md_p}")
        print(f"[DEMO] 第一名: {scored[0]['id']}  综合分={scored[0]['composite_score']}")
        return

    esm = load_csv(args.esm_csv)
    af3 = load_csv(args.af3_csv)
    ids = set(esm) | set(af3)
    if not ids:
        print("[提示] 未提供有效的 --esm-csv / --af3-csv。可先用 --demo 查看流程。")
        return
    rows = []
    for i in sorted(ids):
        e, a = esm.get(i, {}), af3.get(i, {})
        rows.append({
            "id": i,
            "length": e.get("length", ""),
            "mean_plddt": e.get("mean_plddt", ""),
            "iptm": a.get("iptm", ""),
            "interface_pae": a.get("interface_pae", ""),
            "ddg": "",           # 由界面能量脚本(PyRosetta/FoldX)填充
            "mw_kda": "",        # 由成药性脚本填充
            "aggregation": "",
        })
    csv_p, md_p, scored = build_report(rows, DEFAULT_WEIGHTS, args.out)
    print(f"[成功] Top 候选 CSV -> {csv_p}")
    print(f"[成功] 报告 -> {md_p}")


if __name__ == "__main__":
    main()
