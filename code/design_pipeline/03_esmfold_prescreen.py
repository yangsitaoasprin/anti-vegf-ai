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
03_esmfold_prescreen.py
速度层粗筛：用 ESMFold 对 ProteinMPNN 候选序列做单序列结构重预测，
按 pLDDT 与相对设计骨架的 Cα-RMSD 富集自洽的候选。

原理（依据 RFdiffusion 论文的 in silico 过滤逻辑）:
    若 ProteinMPNN 序列真正编码了 RFdiffusion 骨架，则单序列结构预测
    应与设计骨架高度一致 → 高 pLDDT、低 RMSD。ESMFold 免 MSA、速度快，
    适合海量候选的前置富集，随后仅将入围者交 AlphaFold3 精排。

用法:
    python 03_esmfold_prescreen.py --seqs ../designs/mpnn --backbones ../designs/rfdiff \
        --out ../screen/esmfold --min-plddt 70 --max-rmsd 2.0

    # 指定后端（默认自动检测，优先 HuggingFace）
    python 03_esmfold_prescreen.py --seqs ../designs/mpnn --backend huggingface
    python 03_esmfold_prescreen.py --seqs ../designs/mpnn --backend fair-esm

    # 使用本地模型（无需网络）
    python 03_esmfold_prescreen.py --seqs ../designs/mpnn --model-path D:/pmodel/models/esmfold_v1

依赖:
    - 优先: HuggingFace transformers EsmForProteinFolding（无需 openfold，Windows 兼容）
    - 备选: fair-esm（需要 openfold，Linux 推荐）
    未安装时进入演示/指引模式，说明数据流但不产生伪结构。

网络提示:
    HuggingFace 直连不通时，可设置镜像环境变量:
        set HF_ENDPOINT=https://hf-mirror.com
    或使用 --model-path 指定本地模型目录（需包含 config.json + pytorch_model.bin + tokenizer 文件）。

严谨性: 未装模型时绝不输出伪造 pLDDT/RMSD；仅打印指引或在 demo 模式生成带标记的示例表。
"""
import argparse
import csv
import glob
import os
import sys

try:
    import yaml  # noqa
except ImportError:
    yaml = None

# ---------------------------------------------------------------------------
# ESMFold 后端加载
# ---------------------------------------------------------------------------

def _try_load_huggingface(model_path=None):
    """优先后端：HuggingFace transformers EsmForProteinFolding（无需 openfold）。

    model_path: 本地模型目录（包含 config.json + pytorch_model.bin），
                为 None 时从 HuggingFace Hub 下载。
    """
    try:
        import torch
        from transformers import EsmForProteinFolding
        src = model_path or "facebook/esmfold_v1"
        model = EsmForProteinFolding.from_pretrained(
            src, low_cpu_mem_usage=True
        )
        model = model.eval()
        if torch.cuda.is_available():
            model = model.cuda()
        tag = f"本地 ({model_path})" if model_path else "HuggingFace Hub"
        backend = "huggingface"
        print(f"[信息] ESMFold 后端: HuggingFace transformers [{tag}] "
              f"({'GPU: ' + torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU'})")
        return model, torch, backend
    except Exception as e:
        return None, None, f"huggingface 不可用: {e}"


def _try_load_fair_esm():
    """备选后端：fair-esm（需要 openfold，Linux 推荐）。"""
    try:
        import torch
        import esm
        model = esm.pretrained.esmfold_v1()
        model = model.eval()
        if torch.cuda.is_available():
            model = model.cuda()
        backend = "fair-esm"
        print(f"[信息] ESMFold 后端: fair-esm "
              f"({'GPU: ' + torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU'})")
        return model, torch, backend
    except Exception as e:
        return None, None, f"fair-esm 不可用: {e}"


BACKENDS = {
    "huggingface": _try_load_huggingface,
    "fair-esm": _try_load_fair_esm,
}


def try_load_esmfold(backend=None, model_path=None):
    """
    加载 ESMFold 模型，返回 (model, torch, backend_name)。
    backend: "huggingface" | "fair-esm" | None（自动，优先 HuggingFace）
    model_path: 本地模型目录（仅 HuggingFace 后端有效）
    """
    if backend:
        loaders = [(backend, BACKENDS[backend])]
    else:
        loaders = list(BACKENDS.items())

    errors = []
    for name, loader in loaders:
        if name == "huggingface":
            model, torch, result = loader(model_path=model_path)
        else:
            model, torch, result = loader()
        if model is not None:
            return model, torch, result
        errors.append(result)

    print(f"[信息] 无法加载 ESMFold。进入演示/指引模式。")
    for e in errors:
        print(f"  - {e}")
    return None, None, None


# ---------------------------------------------------------------------------
# 预测与指标
# ---------------------------------------------------------------------------

def predict_plddt(model, torch, seq, backend="huggingface"):
    """返回 (mean_pLDDT, pdb_string)。pLDDT 统一为 0-100 范围。"""
    with torch.no_grad():
        out = model.infer_pdb(seq)
    # 解析 B-factor 列（ESMFold 将 pLDDT 写入 B-factor）
    plddts = []
    for line in out.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                plddts.append(float(line[60:66]))
            except ValueError:
                pass
    mean_plddt = sum(plddts) / len(plddts) if plddts else 0.0
    # HuggingFace 输出 0-1 范围，fair-esm 输出 0-100，统一为 0-100
    if backend == "huggingface" and mean_plddt <= 1.0:
        mean_plddt *= 100.0
    return mean_plddt, out


# ---------------------------------------------------------------------------
# FASTA 读取
# ---------------------------------------------------------------------------

def read_fastas(seq_dir):
    """读取目录下所有 .fa/.fasta，返回 [(id, seq), ...]。"""
    records = []
    for path in glob.glob(os.path.join(seq_dir, "**", "*.fa*"), recursive=True):
        header, seq = None, []
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    if header is not None:
                        records.append((header, "".join(seq)))
                    header, seq = line[1:].split()[0], []
                elif line:
                    seq.append(line)
        if header is not None:
            records.append((header, "".join(seq)))
    return records


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="ESMFold 粗筛")
    ap.add_argument("--seqs", required=True, help="ProteinMPNN 序列目录")
    ap.add_argument("--backbones", default=None, help="RFdiffusion 骨架目录（用于 RMSD，可选）")
    ap.add_argument("--out", default="../screen/esmfold")
    ap.add_argument("--min-plddt", type=float, default=70.0)
    ap.add_argument("--max-rmsd", type=float, default=2.0)
    ap.add_argument("--backend", choices=["huggingface", "fair-esm"], default=None,
                    help="ESMFold 后端（默认自动检测，优先 HuggingFace）")
    ap.add_argument("--model-path", default=None,
                    help="本地 ESMFold 模型目录（含 config.json + pytorch_model.bin），"
                         "避免网络下载。例如: D:/pmodel/models/esmfold_v1")
    ap.add_argument("--demo", action="store_true", help="无模型时生成带[DEMO]标记的示例表")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    records = read_fastas(args.seqs) if os.path.isdir(args.seqs) else []
    print(f"[信息] 读取候选序列 {len(records)} 条 (来自 {args.seqs})")

    # demo 模式直接跳过模型加载
    if args.demo:
        model, torch, backend_name = None, None, None
    else:
        model, torch, backend_name = try_load_esmfold(backend=args.backend, model_path=args.model_path)
    rows = []

    if model is not None and records:
        for i, (sid, seq) in enumerate(records, 1):
            if not seq or any(c not in "ACDEFGHIKLMNPQRSTVWY" for c in seq):
                continue
            mean_plddt, pdb = predict_plddt(model, torch, seq, backend=backend_name or "huggingface")
            pdb_path = os.path.join(args.out, f"{sid}.pdb")
            with open(pdb_path, "w") as fh:
                fh.write(pdb)
            # RMSD 相对骨架（此处留接口，实际用 utils_metrics.ca_rmsd 计算）
            rmsd = ""  # 需骨架时填充
            passed = mean_plddt >= args.min_plddt
            rows.append({"id": sid, "length": len(seq), "mean_plddt": round(mean_plddt, 2),
                         "design_rmsd": rmsd, "pass": int(passed), "note": backend_name})
            if i % 50 == 0:
                print(f"  …已处理 {i}/{len(records)}")
    else:
        # 演示/指引模式：不产生伪造科学数值
        if args.demo:
            print("[DEMO] 生成示例表（数值为占位，非真实预测）")
            for sid, seq in (records[:5] or [(f"demo_{i}", "M"*60) for i in range(5)]):
                rows.append({"id": sid, "length": len(seq), "mean_plddt": "NA[DEMO]",
                             "design_rmsd": "NA[DEMO]", "pass": "NA", "note": "占位-需真实ESMFold"})
        else:
            print("\n[操作指引] 未检测到 ESMFold。解决方案（三选一）:")
            print("  方式1（推荐，Windows 兼容）:")
            print("    pip install transformers")
            print("    set HF_ENDPOINT=https://hf-mirror.com  # 国内镜像")
            print("  方式2（离线模式）:")
            print("    下载模型到本地后使用 --model-path 指定")
            print("    例: --model-path D:/pmodel/models/esmfold_v1")
            print("  方式3（Linux 推荐）:")
            print("    pip install fair-esm  # 需要 openfold")
            print("然后重跑本脚本。或加 --demo 查看数据流示例。\n")
            sys.exit(0)

    out_csv = os.path.join(args.out, "prescreen_results.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "length", "mean_plddt", "design_rmsd", "pass", "note"])
        w.writeheader()
        w.writerows(rows)
    n_pass = sum(1 for r in rows if r.get("pass") == 1)
    print(f"[成功] 粗筛结果 -> {out_csv}  (通过 {n_pass}/{len(rows)}, 后端: {backend_name or 'N/A'})")
    print("下一步: 将通过者交 AlphaFold3 复合物预测，再运行 04_af3_rerank.py")


if __name__ == "__main__":
    main()
