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
05b_waltz_aggregation.py
================================================================
用 WALTZ（waltz.switchlab.org，Switch Lab / VIB）对候选微蛋白做淀粉样聚集预测，
输出「风险列」供 05 综合排序使用（对齐 XuL 2026 的 WALTZ 聚集过滤策略）。

WALTZ 简介
----------
WALTZ 是基于 PSSM（位置特异性打分矩阵）+ 物理性质 + 结构伪能量 的淀粉样区域预测器
（Maurer-Stroh et al. 2010, Nat. Methods），输出序列中可能形成 cross-β 淀粉样纤维的
六肽区域。阈值 92 = "Best Overall Performance"（官方推荐默认）。

用法
----
  python 05b_waltz_aggregation.py                          # 默认读 ../results/candidates_6_final.csv
  python 05b_waltz_aggregation.py --in  xx.csv --out 结果.csv

输出
----
  ../results/waltz_aggregation.csv   每个候选 + WALTZ 风险列
  ../results/waltz_aggregation.md    可读报告
"""

import argparse
import csv
import os
import re
import sys
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.normpath(os.path.join(HERE, ".."))
WALTZ_URL = "https://waltz.switchlab.org/results.cgi"
THRESHOLD = "92"          # Best Overall Performance（官方默认）
PH = "7.0"
OUTPUT_MODE = "text_short"  # 结果直接在 HTML 表格中，无需 zip

UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AntiVEGF-miniprotein/1.0"
RETRIES = 4
SLEEP_BASE = 3.0           # 提交间隔/轮询间隔（秒），礼貌访问


def http_post_form(url, fields):
    """application/x-www-form-urlencoded POST，返回 bytes。"""
    body = urllib.parse.urlencode(fields).encode("utf-8")
    req = urllib.request.Request(url, data=body, headers={
        "User-Agent": UA,
        "Content-Type": "application/x-www-form-urlencoded",
    })
    with urllib.request.urlopen(req, timeout=60) as resp:
        return resp.read().decode("utf-8", errors="ignore")


def http_get(url):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=30) as resp:
        return resp.read().decode("utf-8", errors="ignore")


def submit_sequence(name, seq, retries=RETRIES):
    """提交单条 FASTA 到 WALTZ，返回结果 HTML；失败重试。"""
    fields = {
        "sequence": f">{name}\n{seq}",
        "threshold": THRESHOLD,
        "ph": PH,
        "output": OUTPUT_MODE,
        "Submit": "Submit sequences",
    }
    last_err = None
    for attempt in range(retries):
        try:
            html = http_post_form(WALTZ_URL, fields)
            m = re.search(r'href="(OUTPUT/WaltzJob_\d+/WaltzJob_\d+\.html)"', html)
            if not m:
                # 可能服务器忙/校验失败
                if "error" in html.lower() or "correct the errors" in html.lower():
                    raise RuntimeError("WALTZ 服务器拒绝提交（序列含非法字符？）")
                raise RuntimeError("未在响应中找到 job 链接")
            # 轮询等待结果生成
            job_rel = m.group(1)
            job_url = "https://waltz.switchlab.org/" + job_rel
            for _ in range(6):
                time.sleep(SLEEP_BASE)
                try:
                    res = http_get(job_url)
                    if "Regions" in res or "table" in res.lower():
                        return res
                except Exception:
                    pass
            return http_get(job_url)  # 最后再试一次
        except Exception as e:
            last_err = e
            time.sleep(SLEEP_BASE * (attempt + 1))
    raise RuntimeError(f"WALTZ 提交失败（{retries} 次重试后）: {last_err}")


def parse_regions(html, name):
    """从结果 HTML 表格解析该序列的 WALTZ 区域，返回 [(start, end), ...]。"""
    # 提取表格行
    rows = re.findall(r"<tr>(.*?)</tr>", html, re.S)
    for row in rows:
        cells = re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", row, re.S)
        cells = [re.sub(r"<[^>]+>", "", c).strip() for c in cells]
        if len(cells) >= 2 and cells[0] == name:
            reg_str = cells[1].strip()
            if not reg_str:
                return []
            regions = []
            for part in reg_str.split(";"):
                part = part.strip()
                if not part:
                    continue
                mm = re.match(r"^(\d+)-(\d+)$", part)
                if mm:
                    regions.append((int(mm.group(1)), int(mm.group(2))))
            return regions
    return None  # 未找到该序列的行


def risk_grade(regions, seq_len):
    """风险分级：0 区域=低；1-2 短区域=中；>=3 区域或覆盖>10 残基=高。"""
    n = len(regions)
    covered = sum(e - s + 1 for s, e in regions)
    if n == 0:
        return "低", 0
    if n >= 3 or covered >= 11:
        return "高", covered
    return "中", covered


def main():
    ap = argparse.ArgumentParser(description="WALTZ 淀粉样聚集预测 → 风险列")
    ap.add_argument("--in", dest="inp", default=os.path.join(PROJECT, "results", "candidates_6_final.csv"))
    ap.add_argument("--out-csv", default=os.path.join(PROJECT, "results", "waltz_aggregation.csv"))
    ap.add_argument("--out-md", default=os.path.join(PROJECT, "results", "waltz_aggregation.md"))
    args = ap.parse_args()

    if not os.path.exists(args.inp):
        sys.exit(f"[错误] 找不到输入: {args.inp}")

    with open(args.inp, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    seq_field = "binder_seq" if "binder_seq" in rows[0] else "seq"
    print(f"[信息] 读取 {len(rows)} 条候选，序列字段: {seq_field}")

    out_rows = []
    for i, r in enumerate(rows, 1):
        cid = r["id"]
        seq = (r.get(seq_field) or "").strip()
        if not seq:
            print(f"[跳过] {cid} 无序列")
            continue
        print(f"[{i}/{len(rows)}] 提交 {cid} (len={len(seq)}) ...")
        try:
            html = submit_sequence(cid, seq)
        except Exception as e:
            print(f"  [错误] {e}")
            out_rows.append({**r, "waltz_regions": "ERROR", "waltz_n_regions": "ERROR",
                             "waltz_covered_res": "ERROR", "waltz_risk": "ERROR"})
            continue
        regions = parse_regions(html, cid)
        if regions is None:
            print(f"  [警告] 结果页未找到 {cid} 行，按无区域处理")
            regions = []
        grade, covered = risk_grade(regions, len(seq))
        reg_str = "; ".join(f"{s}-{e}" for s, e in regions) or "-"
        print(f"  [结果] 区域={reg_str} | 数量={len(regions)} | 覆盖={covered} aa | 风险={grade}")
        out_rows.append({**r, "waltz_regions": reg_str, "waltz_n_regions": len(regions),
                         "waltz_covered_res": covered, "waltz_risk": grade})
        time.sleep(2)  # 序列间礼貌间隔

    # 写 CSV
    with open(args.out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f"[OK] CSV -> {args.out_csv}")

    # 写 MD
    with open(args.out_md, "w", encoding="utf-8") as f:
        f.write("# WALTZ 淀粉样聚集预测（候选微蛋白）\n\n")
        f.write(f"- 输入: {os.path.basename(args.inp)} | 阈值: Best Overall Performance (92) | pH {PH}\n")
        f.write(f"- WALTZ: Maurer-Stroh et al. 2010 Nat. Methods (waltz.switchlab.org)\n")
        f.write(f"- 风险分级: 低=0 区域; 中=1-2 区域且覆盖≤10aa; 高=≥3 区域或覆盖≥11aa\n\n")
        f.write("| 候选 | 长度 | WALTZ 淀粉样区域 | 区域数 | 覆盖(aa) | 风险 |\n")
        f.write("|---|---|---|---|---|---|\n")
        for r in out_rows:
            f.write(f"| {r['id']} | {r.get('length','')} | {r.get('waltz_regions','-')} | "
                    f"{r.get('waltz_n_regions','-')} | {r.get('waltz_covered_res','-')} | "
                    f"{r.get('waltz_risk','-')} |\n")
        n_high = sum(1 for r in out_rows if r.get("waltz_risk") == "高")
        n_ok = sum(1 for r in out_rows if r.get("waltz_risk") in ("低", "中"))
        f.write(f"\n## 结论\n- 高风险候选: **{n_high}** 个（建议在 05 综合排序中扣除聚集 penalty）\n")
        f.write(f"- 低/中风险候选: **{n_ok}** 个\n")
    print(f"[OK] MD  -> {args.out_md}")


if __name__ == "__main__":
    main()
