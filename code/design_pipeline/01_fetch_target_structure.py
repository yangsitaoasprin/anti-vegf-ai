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
01_fetch_target_structure.py
下载抗 VEGF 课题的靶点结构（VEGF-A / VEGFR 复合物）自 RCSB PDB。

用法:
    python 01_fetch_target_structure.py --pdb 1FLT --out ../data
    python 01_fetch_target_structure.py --pdb 1FLT 3V2A 2FJH --out ../data

说明:
    - 仅使用标准库 urllib，无第三方依赖即可运行（requests 可选）。
    - 下载后打印链、残基数等基本信息，便于后续 hotspot 定义。
    - 严谨性: 脚本不臆造结构；下载失败会明确报错，不生成伪文件。
"""
import argparse
import os
import sys
import urllib.request
import urllib.error

RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb}.pdb"
RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb}.cif"


def download(pdb_id: str, out_dir: str) -> str:
    """下载单个 PDB，优先 .pdb，失败回退 .cif。返回本地文件路径。"""
    pdb_id = pdb_id.strip().upper()
    os.makedirs(out_dir, exist_ok=True)
    for url_tpl, ext in ((RCSB_PDB_URL, "pdb"), (RCSB_CIF_URL, "cif")):
        url = url_tpl.format(pdb=pdb_id)
        dst = os.path.join(out_dir, f"{pdb_id}.{ext}")
        try:
            print(f"[下载] {url}")
            urllib.request.urlretrieve(url, dst)
            print(f"[成功] 已保存 -> {dst}")
            return dst
        except urllib.error.HTTPError as e:
            print(f"[警告] {ext} 格式获取失败 (HTTP {e.code})，尝试下一格式…")
        except Exception as e:  # noqa
            print(f"[警告] 下载 {url} 异常: {e}")
    raise RuntimeError(f"无法下载 {pdb_id}，请检查 PDB 编号或网络。")


def summarize_pdb(path: str) -> None:
    """打印 PDB 文件的链与残基概览（纯文本解析，无需 Biopython）。"""
    if not path.endswith(".pdb"):
        print(f"[信息] {path} 为 CIF 格式，概览请用 Biopython/gemmi 解析。")
        return
    chains = {}
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21]
                resseq = line[22:26].strip()
                chains.setdefault(chain, set()).add(resseq)
    print("[概览] 链 / 残基数：")
    for ch in sorted(chains):
        print(f"    链 {ch}: {len(chains[ch])} 残基")
    print("    ↑ 请据此在 config.yaml 中核对 ligand_chain (VEGF-A) 与 receptor_chain。")


def main():
    ap = argparse.ArgumentParser(description="下载 VEGF 靶点结构 (RCSB PDB)")
    ap.add_argument("--pdb", nargs="+", required=True, help="一个或多个 PDB 编号")
    ap.add_argument("--out", default="../data", help="输出目录")
    args = ap.parse_args()

    ok = 0
    for pid in args.pdb:
        try:
            p = download(pid, args.out)
            summarize_pdb(p)
            ok += 1
        except Exception as e:
            print(f"[错误] {pid}: {e}", file=sys.stderr)
    print(f"\n完成：{ok}/{len(args.pdb)} 个结构下载成功。")


if __name__ == "__main__":
    main()
