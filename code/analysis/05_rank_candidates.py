# =============================================================================
# Author : Yang Sitao (杨四涛)
# Affiliation : College of Pharmacy, Dali University;
#               The Third Affiliated Hospital of Dali University
# Contact : yangsitaoasprin@swpu.edu.cn
# License : MIT (see ../../LICENSE)
# =============================================================================
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05_rank_candidates.py
================================================================
把 O4 产出的能量 CSV（energy_batch.csv）整理成一张「按 ΔG_bind 排序」的候选表。

排序规则
--------
ΔG_bind 越负 = 结合越强 = 排名越靠前（rank 1 最优）。
失败/缺失(ERROR/N/A) 的候选排在最后，不参与最优排名。

输出
----
  results/energy_ranked.csv   带 rank 列的完整排序表
  results/energy_ranked.md    可读的 Markdown 候选表（含 Top10 与统计）

用法
----
  python 05_rank_candidates.py --in  ../results/energy_batch.csv \
      --out-csv ../results/energy_ranked.csv --out-md ../results/energy_ranked.md
"""
import argparse
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.normpath(os.path.join(HERE, ".."))


def parse_float(v):
    """解析数值；失败返回 None。"""
    if v is None:
        return None
    s = str(v).strip().replace("kcal/mol", "").strip()
    if s in ("", "N/A", "ERROR", "None"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


def main():
    ap = argparse.ArgumentParser(description="O4 能量结果 → ΔG 排序候选表")
    ap.add_argument("--in", dest="inp", default=os.path.join(PROJECT, "results", "energy_batch.csv"))
    ap.add_argument("--out-csv", default=os.path.join(PROJECT, "results", "energy_ranked.csv"))
    ap.add_argument("--out-md", default=os.path.join(PROJECT, "results", "energy_ranked.md"))
    ap.add_argument("--top", type=int, default=10, help="Markdown 表中展示的 TopN")
    args = ap.parse_args()

    if not os.path.exists(args.inp):
        print(f"[错误] 找不到输入: {args.inp}")
        raise SystemExit(1)

    rows = []
    with open(args.inp, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            rows.append(r)

    # 计算数值 ΔG 用于排序
    for r in rows:
        r["_dg"] = parse_float(r.get("prime_dg"))
        r["_vdw"] = parse_float(r.get("dg_vdw"))
        r["_hbond"] = parse_float(r.get("dg_hbond"))
        r["_solv"] = parse_float(r.get("dg_solv_gb"))

    # 排序：有有效 ΔG 的在前（升序，最负最优），无效(无)的排最后
    def sort_key(r):
        dg = r["_dg"]
        if dg is None:
            return (1, 0.0)        # 无效组，靠后
        return (0, dg)             # 有效组，按 ΔG 升序

    rows.sort(key=sort_key)

    # 赋 rank（仅有效候选连续编号；无效候选 rank 标 "-")
    fields = ["id", "length", "prime_dg", "dg_vdw", "dg_hbond", "dg_solv_gb",
              "mean_plddt", "docked", "pass", "note"]
    out_rows = []
    rank = 0
    for r in rows:
        if r["_dg"] is not None:
            rank += 1
            r["rank"] = str(rank)
        else:
            r["rank"] = "-"
        out_rows.append(r)

    # 写 CSV（带 rank 列）
    csv_fields = ["rank"] + fields
    with open(args.out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=csv_fields, extrasaction="ignore")
        w.writeheader()
        for r in out_rows:
            w.writerow(r)
    print(f"[OK] 排序表 CSV -> {args.out_csv}  ({len(out_rows)} 行)")

    # 统计
    valid = [r for r in out_rows if r["_dg"] is not None]
    best = valid[0] if valid else None
    worst = valid[-1] if valid else None
    n_valid = len(valid)
    n_fail = len(out_rows) - n_valid

    # 写 Markdown
    with open(args.out_md, "w", encoding="utf-8") as f:
        f.write("# 抗 VEGF-A 微蛋白 · ΔG 排序候选表\n\n")
        f.write(f"- 候选总数: **{len(out_rows)}** ｜ 有效评估: **{n_valid}** ｜ 计算失败: **{n_fail}**\n")
        if best:
            f.write(f"- **最优候选 (rank 1)**: `{best['id']}` — ΔG_bind = **{best['prime_dg']} kcal/mol** "
                    f"(vdW {best.get('dg_vdw','')} / Hbond {best.get('dg_hbond','')} / Solv_GB {best.get('dg_solv_gb','')})\n")
        f.write(f"- 排序规则: ΔG_bind 越负结合越强，rank 1 最优；ERROR/N/A 排末尾。\n\n")

        f.write(f"## Top {min(args.top, n_valid)} 候选\n\n")
        f.write("| Rank | ID | Len | ΔG_bind (kcal/mol) | vdW | Hbond | Solv_GB | 备注 |\n")
        f.write("|-----:|----|----:|--------------------:|----:|------:|--------:|------|\n")
        for r in valid[:args.top]:
            f.write(f"| {r['rank']} | `{r['id']}` | {r.get('length','')} | **{r['prime_dg']}** | "
                    f"{r.get('dg_vdw','')} | {r.get('dg_hbond','')} | {r.get('dg_solv_gb','')} | {r.get('note','')} |\n")

        if n_fail:
            f.write(f"\n## 计算失败 / 无效 ({n_fail})\n\n")
            f.write("| ID | 备注 |\n|----|------|\n")
            for r in out_rows:
                if r["_dg"] is None:
                    f.write(f"| `{r['id']}` | {r.get('note','')} |\n")

    print(f"[OK] 排序表 MD  -> {args.out_md}")


if __name__ == "__main__":
    main()
