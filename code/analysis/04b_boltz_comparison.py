# =============================================================================
# Author : Yang Sitao (杨四涛)
# Affiliation : School of Pharmaceutical Sciences, Dali University;
#               The Third Affiliated Hospital of Dali University
# Contact : yangsitaoasprin@126.com
# License : MIT (see ../../LICENSE)
# =============================================================================
#!/usr/bin/env python3
"""
04b_boltz_comparison.py
Boltz-1 模型对照分析：比较 AF3 与 Boltz-1 对相同候选的复合物预测一致性。
用于验证双层筛选策略的模型独立性。

用法：
    python 04b_boltz_comparison.py --af3-csv ../results/af3_rerank.csv --boltz-csv ../results/boltz_results.csv
"""
import argparse, csv, os, math

def spearman_rho(xs, ys):
    n = len(xs)
    if n < 3: return 0.0, 1.0
    rx = [sorted(xs).index(v) for v in xs]
    ry = [sorted(ys).index(v) for v in ys]
    d = sum((rx[i] - ry[i])**2 for i in range(n))
    rho = 1 - 6*d/(n*(n*n-1))
    # 近似 p 值（无需 scipy）
    t = rho * math.sqrt((n-2)/(1-rho*rho)) if abs(rho)<1 else float('inf')
    import math as _math
    if t == float('inf'):
        p = 0.0
    else:
        # 用正态分布近似
        p = _math.erfc(abs(t)/_math.sqrt(2))
    return rho, min(p, 1.0)

def kendall_tau(r1, r2):
    """排名向量的一致性"""
    n = len(r1)
    if n < 2: return 0.0, 1.0
    conc, disc = 0, 0
    for i in range(n):
        for j in range(i+1, n):
            if (r1[i]-r1[j])*(r2[i]-r2[j]) > 0:
                conc += 1
            else:
                disc += 1
    tau = (conc - disc) / (conc + disc) if (conc + disc) > 0 else 0.0
    var = (4*n+10)/(9*n*(n-1)) if n > 1 else 1
    z = tau / math.sqrt(var) if var > 0 else 0
    p = math.erfc(abs(z)/math.sqrt(2)) if z != 0 else 1.0
    return tau, min(p, 1.0)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--af3-csv", help="AF3 精排结果 CSV")
    p.add_argument("--boltz-csv", help="Boltz-1 结果 CSV")
    p.add_argument("--out", default="../results/boltz_comparison_report.md")
    p.add_argument("--demo", action="store_true", help="生成模拟数据验证流程")
    args = p.parse_args()

    if args.demo:
        print("[demo] Boltz-1 vs AF3 一致分析流程验证")
        print("[demo] 模拟 20 个候选的 ipTM 数据...")
        import random
        random.seed(42)
        ids = [f"cand_{i:03d}" for i in range(20)]
        xs = [random.uniform(0.3, 0.9) for _ in range(20)]
        ys = [x * random.gauss(1.0, 0.1) for x in xs]  # 高相关模拟
        ys = [min(0.95, max(0.2, y)) for y in ys]
        rho, _ = spearman_rho(xs, ys)
        r1 = sorted(range(len(xs)), key=lambda i: xs[i])
        r2 = sorted(range(len(ys)), key=lambda i: ys[i])
        tau, _ = kendall_tau(r1, r2)
        print(f"  模拟 Spearman rho = {rho:.4f}")
        print(f"  模拟 Kendall tau = {tau:.4f}")
        print(f"  判断: {'模型独立' if tau>=0.6 else '模型敏感'}")
        print("[demo] 完成。准备真实数据后运行: python 04b_boltz_comparison.py --af3-csv ... --boltz-csv ...")
        return

    args = p.parse_args()

    # 读取 AF3 结果
    af3 = {}
    with open(args.af3_csv) as f:
        for row in csv.DictReader(f):
            af3[row["id"]] = {"ipTM": float(row.get("ipTM",0)), "pAE": float(row.get("interface_pAE",10))}

    # 读取 Boltz-1 结果
    boltz = {}
    with open(args.boltz_csv) as f:
        for row in csv.DictReader(f):
            boltz[row["id"]] = {"ipTM": float(row.get("ipTM",0)), "pAE": float(row.get("interface_pAE",10))}

    # 取交集
    common = set(af3.keys()) & set(boltz.keys())
    print(f"AF3 候选: {len(af3)}, Boltz-1: {len(boltz)}, 交集: {len(common)}")

    xs = [af3[i]["ipTM"] for i in common]
    ys = [boltz[i]["ipTM"] for i in common]

    rho, p_rho = spearman_rho(xs, ys)
    # 排名
    r1 = [sorted(xs).index(v) for v in xs]
    r2 = [sorted(ys).index(v) for v in ys]
    tau, p_tau = kendall_tau(r1, r2)

    report = [
        "# Boltz-1 vs AF3 精排一致性分析报告",
        "",
        f"比较候选数: {len(common)}",
        f"",
        f"## ipTM 相关性",
        f"- Spearman rho = {rho:.4f} (p = {p_rho:.6f})",
        f"- Kendall tau = {tau:.4f} (p = {p_tau:.6f})",
        f"",
        f"## 判断",
        f"- tau >= 0.6 -> 模型独立性成立：双层筛选结论不依赖于精排模型",
        f"- tau < 0.6 -> 筛选结论对模型敏感，需进一步分析差异原因",
    ]
    verdict = "模型独立" if tau >= 0.6 else "模型敏感"
    report.append(f"")
    report.append(f"## 结论：{verdict}")

    out_path = args.out
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        f.write("\n".join(report))
    print(f"报告已保存: {out_path}")

    # 打印 Top 10 对比
    top10 = sorted(common, key=lambda i: af3[i]["ipTM"], reverse=True)[:10]
    print("\nTop 10 (AF3 排序):")
    print(f"{'ID':20s} {'AF3 ipTM':10s} {'Boltz ipTM':10s} {'AF3 pAE':10s} {'Boltz pAE':10s}")
    for i in top10:
        print(f"{i:20s} {af3[i]['ipTM']:.4f}     {boltz[i]['ipTM']:.4f}     {af3[i]['pAE']:.2f}      {boltz[i]['pAE']:.2f}")

if __name__ == "__main__":
    main()
