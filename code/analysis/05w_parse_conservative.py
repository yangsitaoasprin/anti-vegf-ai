# =============================================================================
# Author : Yang Sitao (杨四涛)
# Affiliation : College of Pharmacy, Dali University;
#               The Third Affiliated Hospital of Dali University
# Contact : yangsitaoasprin@swpu.edu.cn
# License : MIT (see ../../LICENSE)
# =============================================================================
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""05w — 解析保守 #2 候选的 Boltz confidence/pae（避开 glob 对方括号的字符类误判）。
直接 walk 文件系统定位 confidence_*.json / pae_*.npz，复用 05k 的数值解析与界面pAE。
"""
import os, csv, json, sys, importlib.util
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, ".."))
SPEC = importlib.util.spec_from_file_location(
    "bk", os.path.join(HERE, "05k_boltz_redesign_2.py"))
bk = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bk)

IN_FA = os.path.join(PROJ, "designs", "candidates_6", "deA_redesign_1_2",
                     "conservative_if_2", "esmfold_passed.fa")
WORK_ROOT = bk.WORK_ROOT  # results/deA_redesign_1_2/boltz_pdb
SEED = 0
OUT_CSV = os.path.join(PROJ, "designs", "candidates_6", "deA_redesign_1_2",
                       "conservative_if_2", "boltz_s0_parsed.csv")


def walk_find(out_dir, want_base):
    """want_base 形如 confidence_<stem>_model_0 ; 返回完整路径，避免 glob 方括号问题。"""
    hits = []
    for root, _, files in os.walk(out_dir):
        for f in files:
            name = f if isinstance(f, str) else f.decode("utf-8", "replace")
            if name.startswith(want_base) and (
                    name.endswith(".json") or name.endswith(".npz")):
                hits.append(os.path.join(root, f))
    return sorted(hits)


def main():
    seqs = bk.read_fastas(IN_FA)
    rows = []
    for si, (name, seq) in enumerate(seqs, 1):
        stem = f"{name}_s{SEED}"
        out_dir = os.path.join(WORK_ROOT, stem, "out")
        cj = walk_find(out_dir, f"confidence_{stem}_model_0")
        pz = walk_find(out_dir, f"pae_{stem}_model_0")
        if not cj:
            print(f"[{si}/{len(seqs)}] {name}: 无 confidence", flush=True)
            rows.append({"id": name, "len": len(seq), "iptm_mean": "NA",
                         "iptm_best": "NA", "ptm_mean": "NA",
                         "interface_pae_mean": "NA", "interface_plddt_mean": "NA",
                         "n_seed": 1, "pass": 0, "tag": "conservative_if_s0"})
            continue
        d = json.load(open(cj[0]))
        iptm = d.get("iptm"); ptm = d.get("ptm")
        iplddt = d.get("complex_iplddt") if d.get("complex_iplddt") is not None \
            else d.get("complex_plddt")
        ipae = bk.interface_pae_from_npz(pz[0], bk.LEN_A) if pz else "NA"
        iplddt = bk.norm100(iplddt) if iplddt is not None else None
        ok = int(iptm >= 0.5) if isinstance(iptm, (int, float)) else 0
        print(f"[{si}/{len(seqs)}] {name}: ipTM={iptm} ptm={ptm} "
              f"pAE={ipae} pLDDT={iplddt}", flush=True)
        rows.append({"id": name, "len": len(seq),
                     "iptm_mean": iptm if iptm is not None else "NA",
                     "iptm_best": iptm if iptm is not None else "NA",
                     "ptm_mean": ptm if ptm is not None else "NA",
                     "interface_pae_mean": ipae if ipae is not None else "NA",
                     "interface_plddt_mean": iplddt if iplddt is not None else "NA",
                     "n_seed": 1, "pass": ok, "tag": "conservative_if_s0"})

    def sk(r):
        ip = r["iptm_mean"]; pa = r["interface_pae_mean"]
        ip = ip if isinstance(ip, (int, float)) else -1
        pa = pa if isinstance(pa, (int, float)) else 999
        return (-ip, pa)
    rows.sort(key=sk)
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["id", "len", "iptm_mean", "iptm_best",
                              "ptm_mean", "interface_pae_mean", "interface_plddt_mean",
                              "n_seed", "pass", "tag"])
        w.writeheader(); w.writerows(rows)
    print(f"\n[成功] -> {OUT_CSV}")
    print("=== 排名(ipTM降序 / 界面pAE升序) ===")
    for r in rows:
        print(f"  {str(r['id'])[:40]:40s} ipTM={r['iptm_mean']} ptm={r['ptm_mean']} "
              f"pAE={r['interface_pae_mean']} pLDDT={r['interface_plddt_mean']} pass={r['pass']}")
    print(f"\n基线 #2 最优: ipTM={bk.BASE_IPTM} 界面pAE={bk.BASE_IPAE} 界面pLDDT={bk.BASE_IPLDDT}")


if __name__ == "__main__":
    main()
