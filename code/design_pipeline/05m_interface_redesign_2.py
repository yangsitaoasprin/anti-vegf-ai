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
05m_interface_redesign_2.py — #2_T0.2_s83 界面定向进化
以 s83 的 Boltz 好模式复合物(s42)为模板：
  - 靶标链 A (98aa) 全固定
  - 结合剂链 B (100aa) 仅重设计 <=5A 的 31 个界面接触残基，其余(69个)锁定为 s83 身份
  - 该 scropfold 折叠保持不变，只进化结合面 -> 争取 ipTM 0.6+
ProteinMPNN 多温度采样 -> WALTZ 聚集过滤 -> 输出。
"""
import os, sys, subprocess, json, glob, re, shutil, tempfile
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, ".."))
MPNN_DIR = r"D:\ProteinMPNN"
ESM_PY = r"C:\Users\Administrator\.conda\envs\esm\python.exe"

TEMPLATE_PDB = os.path.join(PROJ, "results", "deA_redesign_1_2", "boltz_pdb",
                            "#2_T0.2_s83_s42_model_0.pdb")  # chain A=target, B=binder
TARGET_LEN = 98
BINDER_LEN = 100
PARENT_SEQ = ("LEALKEQLLKYAKTAEKCAKEAQAKAEEEQKKSEYYKKIAEDEKLAEEFAKEIGGISVEE"
              "AKEIAKDLAEAYKLAAEHYKRTSEACQARAEALKKQAAAL")  # s83 设计序列
IFACE_CUT = 5.0  # 界面定义: 残基最小重原子距离 <= 5.0A
TEMPS = [0.1, 0.15, 0.2, 0.3]
N_SEQ = 30
OUT_ROOT = os.path.join(PROJ, "designs", "candidates_6", "deA_redesign_1_2", "interface_evo_2")
os.makedirs(OUT_ROOT, exist_ok=True)

def build_env():
    envdir = os.path.dirname(ESM_PY)
    lib_bin = os.path.join(envdir, "Library", "bin")
    scripts = os.path.join(envdir, "Scripts")
    mingw = os.path.join(envdir, "Library", "mingw-w64", "bin")
    sys32, sysroot = r"C:\Windows\System32", r"C:\Windows"
    pp = [p for p in [envdir, scripts, lib_bin, mingw, sys32, sysroot] if os.path.isdir(p)]
    return {"PATH": ";".join(pp), "SYSTEMROOT": sysroot, "WINDIR": sysroot,
            "COMSPEC": os.path.join(sys32, "cmd.exe"),
            "TEMP": os.environ.get("TEMP", r"C:\Windows\Temp"),
            "USERNAME": os.environ.get("USERNAME", "Administrator"),
            "OS": "Windows_NT", "PROCESSOR_ARCHITECTURE": "AMD64",
            "NUMBER_OF_PROCESSORS": os.environ.get("NUMBER_OF_PROCESSORS", "8"),
            "APPDATA": os.environ.get("APPDATA", ""), "PYTHONIOENCODING": "utf-8"}

def run_cmd(cmd, env, cwd=None, timeout=1200):
    print("  $ " + " ".join(str(c) for c in cmd)[:140], flush=True)
    return subprocess.run(cmd, env=env, cwd=cwd, capture_output=True, text=True, timeout=timeout)

def compute_interface(pdb):
    p = PDBParser(QUIET=True).get_structure("x", pdb)
    A, B = p[0]["A"], p[0]["B"]
    tgt = [a for a in A.get_atoms() if a.element != "H"]
    ns = NeighborSearch(tgt)
    iface = []
    for i, r in enumerate(B):
        rs = [a for a in r if a.element != "H"]
        if not rs:
            continue
        best = 1e9
        for a in rs:
            for o in ns.search(a.coord, IFACE_CUT, level="A"):
                d = float(np.linalg.norm(a.coord - o.coord))
                if d < best:
                    best = d
        if best < 1e9:
            iface.append(i + 1)
    return sorted(iface)

def simple_agg_score(seq):
    a6 = sum(1 for i in range(len(seq)-5) if seq[i:i+6].count("A") >= 4)
    consec = max((len(g) for g in "".join("A" if c=="A" else " " for c in seq).split()), default=0)
    hydro = sum(1 for i in range(len(seq)-5) if sum(seq[i:i+6].count(aa) for aa in "AILMFWV") >= 5)
    return a6 * 10 + consec * 5 + hydro * 2

def run_waltz(seq):
    import urllib.request, urllib.parse, time
    WALTZ_URL = "https://waltz.switchlab.org/results.cgi"
    UA = "Mozilla/5.0 AntiVEGF-miniprotein/1.0"
    fields = {"sequence": f">test\n{seq}", "threshold": "92", "ph": "7.0", "output": "text_short", "Submit": "Submit sequences"}
    try:
        body = urllib.parse.urlencode(fields).encode("utf-8")
        req = urllib.request.Request(WALTZ_URL, data=body, headers={"User-Agent": UA, "Content-Type": "application/x-www-form-urlencoded"})
        with urllib.request.urlopen(req, timeout=60) as resp:
            html = resp.read().decode("utf-8", errors="ignore")
        m = re.search(r'href="(OUTPUT/WaltzJob_\d+/WaltzJob_\d+\.html)"', html)
        if not m:
            return None
        job_url = "https://waltz.switchlab.org/" + m.group(1)
        for _ in range(6):
            time.sleep(2)
            try:
                req2 = urllib.request.Request(job_url, headers={"User-Agent": UA})
                with urllib.request.urlopen(req2, timeout=30) as resp2:
                    res = resp2.read().decode("utf-8", errors="ignore")
                if "Regions" in res:
                    rows = re.findall(r"<tr>(.*?)</tr>", res, re.S)
                    for row in rows:
                        cells = re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", row, re.S)
                        cells = [re.sub(r"<[^>]+>", "", c).strip() for c in cells]
                        if len(cells) >= 2 and cells[0] == "test":
                            return cells[1].strip() == ""
            except Exception:
                pass
        return None
    except Exception:
        return None

def main():
    work = tempfile.mkdtemp(prefix="mpnn_iface_")
    pdb_dir = os.path.join(work, "pdbs"); os.makedirs(pdb_dir)
    pdb_name = os.path.splitext(os.path.basename(TEMPLATE_PDB))[0]
    shutil.copy(TEMPLATE_PDB, os.path.join(pdb_dir, pdb_name + ".pdb"))

    iface = compute_interface(TEMPLATE_PDB)
    print(f"[Interface] <= {IFACE_CUT}A contacting residues (1-based): {iface} (n={len(iface)})")
    fixed_B = [i for i in range(1, BINDER_LEN + 1) if i not in iface]
    print(f"[Fix] B chain fixed (non-interface) = {len(fixed_B)} pos; designable = {len(iface)}")

    fix_path = os.path.join(work, "fixed_positions.jsonl")
    # 占位, parse 后改写为真实 parsed name
    with open(fix_path, "w") as f:
        json.dump({pdb_name: {"A": list(range(1, TARGET_LEN + 1)), "B": fixed_B}}, f); f.write("\n")

    env = build_env()
    parse_py = os.path.join(MPNN_DIR, "helper_scripts", "parse_multiple_chains.py")
    assign_py = os.path.join(MPNN_DIR, "helper_scripts", "assign_fixed_chains.py")
    run_py = os.path.join(MPNN_DIR, "protein_mpnn_run.py")
    weights = os.path.join(MPNN_DIR, "vanilla_model_weights")

    print("[Step1] parse chains")
    parsed = os.path.join(work, "parsed_chains.jsonl")
    r = run_cmd([ESM_PY, parse_py, f"--input_path={pdb_dir}", f"--output_path={parsed}"], env)
    if r.returncode != 0:
        print("FAIL parse:", r.stderr[-300:]); return
    # 修正 parsed name -> 重写 fix_path
    if os.path.exists(parsed):
        with open(parsed) as f:
            raw_name = json.loads(f.readline()).get("name", pdb_name)
        with open(fix_path, "w") as f:
            json.dump({raw_name: {"A": list(range(1, TARGET_LEN + 1)), "B": fixed_B}}, f); f.write("\n")
        print(f"  parsed name = {raw_name}")

    print("[Step2] assign chains (design B only)")
    chain_id = os.path.join(work, "chain_id_map.jsonl")
    r = run_cmd([ESM_PY, assign_py, f"--input_path={parsed}", f"--output_path={chain_id}", "--chain_list=B"], env)
    if r.returncode != 0:
        print("FAIL assign:", r.stderr[-300:]); return

    all_seqs = {}
    for temp in TEMPS:
        t_out = os.path.join(work, f"mpnn_out_{temp}")
        os.makedirs(t_out, exist_ok=True)
        print(f"[Step3] MPNN T={temp}")
        r = run_cmd([ESM_PY, run_py, "--jsonl_path", parsed, "--chain_id_jsonl", chain_id,
                     "--fixed_positions_jsonl", fix_path, "--out_folder", t_out,
                     "--num_seq_per_target", str(N_SEQ), "--sampling_temp", str(temp),
                     "--batch_size", "1", "--seed", "42",
                     "--path_to_model_weights", weights, "--model_name", "v_48_020"], env, timeout=1800)
        if r.returncode != 0:
            print(f"  FAIL T={temp}: {r.stderr[-300:]}"); continue
        seq_dir = os.path.join(t_out, "seqs")
        for fa in glob.glob(os.path.join(seq_dir, "*.fa")):
            with open(fa) as f:
                lines = [l.strip() for l in f.readlines()]
            for i in range(2, len(lines)-1, 2):  # 跳过第一组(A链=靶标)
                if lines[i].startswith(">") and lines[i+1]:
                    seq = lines[i+1]
                    if len(seq) == BINDER_LEN:
                        name = f"#2if_T{temp}_s{len(all_seqs)}"
                        all_seqs[name] = seq
        print(f"  T={temp}: accumulated {len(all_seqs)} seqs")

    if not all_seqs:
        print("No sequences generated"); return

    # 预筛: A-rich 密度
    print(f"[Step4] pre-filter {len(all_seqs)} seqs (agg score, keep top 40)")
    scored = []
    for name, seq in all_seqs.items():
        if len(seq) != BINDER_LEN:
            continue
        scored.append((name, seq, simple_agg_score(seq)))
    scored.sort(key=lambda x: x[2])
    top = scored[:40]
    print(f"  agg range {top[0][2]}-{top[-1][2]}")

    # WALTZ
    print(f"[Step5] WALTZ check top {len(top)}")
    passed = []
    for i, (name, seq, score) in enumerate(top):
        print(f"  [{i+1}/{len(top)}] {name}", end=" ", flush=True)
        ok = run_waltz(seq)
        if ok:
            passed.append((name, seq, score)); print("PASS")
        elif ok is False:
            print("FAIL")
        else:
            print("ERR")
    print(f"[Done] WALTZ passed: {len(passed)}/{len(top)}")

    out_fa = os.path.join(OUT_ROOT, "interface_evo_WALTZ_passed.fa")
    with open(out_fa, "w") as f:
        for name, seq, score in passed:
            muts = sum(1 for a, b in zip(PARENT_SEQ, seq) if a != b)
            a_new = seq.count("A")
            f.write(f">{name} muts={muts} A={a_new} agg_score={score}\n{seq}\n")
    print(f"OUTPUT: {out_fa}")

if __name__ == "__main__":
    main()
