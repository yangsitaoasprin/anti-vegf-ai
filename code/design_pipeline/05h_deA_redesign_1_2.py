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
05h_deA_redesign_1_2.py — #1 和 #2 的去 A 重设计
复用 05d 的 ProteinMPNN 调用链（PDB解析 + assign_fixed_chains + WALTZ）。
加速: ProteinMPNN → 简单聚集评分预筛(top 30) → WALTZ → 输出。

用法: python 05h_deA_redesign_1_2.py
"""
import os, sys, subprocess, json, glob, re, shutil, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.normpath(os.path.join(HERE, ".."))
MPNN_DIR = r"D:\ProteinMPNN"
ESM_PY = r"C:\Users\Administrator\.conda\envs\esm\python.exe"

# ============ 候选配置 ============
CONFIGS = {
    "#1": {
        "pdb": os.path.join(PROJ, "designs", "candidates_6", "structures", "vegf_len80_1_T0.2_s2.complex.pdb"),
        "binder_chain": "B",
        "seq_orig": "SPGLELCRRAVERLRACADAIAAAGTLEAARACAASLPADSAPLAEALVANAAAGTAAETAAALRAAAEALERLCAEQEA",
        "fixed_pos": [2,3,6,7,10,14,40,43,44,47,48,51,53,54,55,59,60,63,64,67,70,71,72],
    },
    "#2": {
        "pdb": os.path.join(PROJ, "designs", "candidates_6", "structures", "vegf_len100_4_T0.2_s1.complex.pdb"),
        "binder_chain": "B",
        "seq_orig": "KEAEKAALLAQAAAAQAAAATAQAKAAEAQAKAEEYKAIAENEDLAKEFAKKIGNISVEEAKEIAKDLAEAYRLAAEHHQQTAAACQARADALLAQAKAL",
        "fixed_pos": [53,54,56,60,63,64,67,68,70,71,74,75,78,82,85,86,88,89,92,93,96],
    },
}
TEMPS = [0.1, 0.2, 0.3]
N_SEQ = 50
OUT_ROOT = os.path.join(PROJ, "designs", "candidates_6", "deA_redesign_1_2")
os.makedirs(OUT_ROOT, exist_ok=True)

def build_env():
    envdir = os.path.dirname(ESM_PY)
    lib_bin = os.path.join(envdir, "Library", "bin")
    scripts = os.path.join(envdir, "Scripts")
    mingw = os.path.join(envdir, "Library", "mingw-w64", "bin")
    sys32, sysroot = r"C:\Windows\System32", r"C:\Windows"
    pp = [p for p in [envdir, scripts, lib_bin, mingw, sys32, sysroot] if os.path.isdir(p)]
    return {
        "PATH": ";".join(pp), "SYSTEMROOT": sysroot, "WINDIR": sysroot,
        "COMSPEC": os.path.join(sys32, "cmd.exe"),
        "TEMP": os.environ.get("TEMP", r"C:\Windows\Temp"),
        "USERNAME": os.environ.get("USERNAME", "Administrator"),
        "OS": "Windows_NT", "PROCESSOR_ARCHITECTURE": "AMD64",
        "NUMBER_OF_PROCESSORS": os.environ.get("NUMBER_OF_PROCESSORS", "8"),
        "APPDATA": os.environ.get("APPDATA", ""),
        "PYTHONIOENCODING": "utf-8",
    }

def run_cmd(cmd, env, cwd=None, timeout=600):
    print("  $ " + " ".join(str(c) for c in cmd)[:120], flush=True)
    return subprocess.run(cmd, env=env, cwd=cwd, capture_output=True, text=True, timeout=timeout)

def simple_agg_score(seq):
    """快速聚集评分: A-rich 6mer密度 + 连续A长度"""
    a6 = sum(1 for i in range(len(seq)-5) if seq[i:i+6].count("A") >= 4)
    consec = max((len(g) for g in "".join("A" if c=="A" else " " for c in seq).split()), default=0)
    hydro = sum(1 for i in range(len(seq)-5) if sum(seq[i:i+6].count(aa) for aa in "AILMFWV") >= 5)
    return a6 * 10 + consec * 5 + hydro * 2

def run_waltz(seq):
    """调用 05b_waltz_aggregation 的 WALTZ API"""
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
        if not m: return None
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
                            return cells[1].strip() == ""  # True = 无区域 (通过)
            except: pass
        return None
    except: return None

def main():
    for cand_id, cfg in CONFIGS.items():
        print(f"\n{'='*60}\n  {cand_id}: fixed {len(cfg['fixed_pos'])} residues\n{'='*60}")
        
        work = tempfile.mkdtemp(prefix=f"mpnn_{cand_id}_")
        pdb_dir = os.path.join(work, "pdbs"); os.makedirs(pdb_dir)
        pdb_name = os.path.splitext(os.path.basename(cfg["pdb"]))[0]
        shutil.copy(cfg["pdb"], os.path.join(pdb_dir, pdb_name + ".pdb"))
        
        # fixed_positions.jsonl — A链固定全部残基(不设计), B链固定界面残基
        fix_path = os.path.join(work, "fixed_positions.jsonl")
        # A链: 全部残基固定 (1-N), B链: 只固定界面残基
        with open(fix_path, "w") as f:
            json.dump({pdb_name: {"A": list(range(1, 99)), cfg["binder_chain"]: cfg["fixed_pos"]}}, f)
            f.write("\n")
        
        env = build_env()
        parse_py = os.path.join(MPNN_DIR, "helper_scripts", "parse_multiple_chains.py")
        assign_py = os.path.join(MPNN_DIR, "helper_scripts", "assign_fixed_chains.py")
        run_py = os.path.join(MPNN_DIR, "protein_mpnn_run.py")
        weights = os.path.join(MPNN_DIR, "vanilla_model_weights")
        
        # Parse chains
        parsed = os.path.join(work, "parsed_chains.jsonl")
        print("[Step 1] parse chains")
        r = run_cmd([ESM_PY, parse_py, f"--input_path={pdb_dir}", f"--output_path={parsed}"], env)
        if r.returncode != 0:
            print(f"  FAIL: {r.stderr[-200:]}")
            continue
        
        # Fix parsed name (Windows rfind path bug)
        if os.path.exists(parsed):
            with open(parsed) as f:
                raw_name = json.loads(f.readline()).get("name", pdb_name)
            with open(fix_path, "w") as f:
                json.dump({raw_name: {"A": list(range(1, 99)), cfg["binder_chain"]: cfg["fixed_pos"]}}, f)
                f.write("\n")
        
        # Assign chains — 只设计 B 链
        chain_id = os.path.join(work, "chain_id_map.jsonl")
        print("[Step 2] assign chains (design B chain only)")
        r = run_cmd([ESM_PY, assign_py, f"--input_path={parsed}", f"--output_path={chain_id}",
                     f"--chain_list={cfg['binder_chain']}"], env)
        if r.returncode != 0:
            print(f"  FAIL: {r.stderr[-200:]}")
            continue
        
        # Run ProteinMPNN (3 temps)
        all_seqs = {}
        for temp in TEMPS:
            t_out = os.path.join(work, f"mpnn_out_{temp}")
            os.makedirs(t_out, exist_ok=True)
            print(f"[Step 3] MPNN T={temp}")
            r = run_cmd([ESM_PY, run_py,
                "--jsonl_path", parsed, "--chain_id_jsonl", chain_id,
                "--fixed_positions_jsonl", fix_path, "--out_folder", t_out,
                "--num_seq_per_target", str(N_SEQ), "--sampling_temp", str(temp),
                "--batch_size", "1", "--seed", "42",
                "--path_to_model_weights", weights, "--model_name", "v_48_020"], env, timeout=1200)
            if r.returncode != 0:
                print(f"  FAIL T={temp}: {r.stderr[-200:]}")
                continue
            # Read output — 多序列 FASTA (A链VEGF+B链binder), 只取B链
            seq_dir = os.path.join(t_out, "seqs")
            for fa in glob.glob(os.path.join(seq_dir, "*.fa")):
                with open(fa) as f:
                    lines = [l.strip() for l in f.readlines()]
                # 每2行为一组 (header, seq), 第一组是A链(VEGF固定), 之后是B链(设计)
                for i in range(2, len(lines)-1, 2):  # 跳过前2行(A链)
                    if lines[i].startswith(">") and lines[i+1]:
                        seq = lines[i+1]
                        if len(seq) == len(cfg["seq_orig"]):
                            name = f"{cand_id}_T{temp}_s{len(all_seqs)}"
                            all_seqs[name] = seq
            print(f"  T={temp}: total seqs={len(all_seqs)}")
        
        if not all_seqs:
            print(f"  No sequences generated for {cand_id}")
            continue
        
        # Pre-filter: simple agg score (no WSL needed)
        print(f"[Step 4] Pre-filter {len(all_seqs)} seqs by aggregation score")
        scored = []
        for name, seq in all_seqs.items():
            if len(seq) != len(cfg["seq_orig"]): continue
            s = simple_agg_score(seq)
            scored.append((name, seq, s))
        scored.sort(key=lambda x: x[2])
        top30 = scored[:30]
        print(f"  Kept top 30 (score range {top30[0][2]}-{top30[-1][2]})")
        
        # WALTZ check
        print(f"[Step 5] WALTZ check top 30")
        passed = []
        for i, (name, seq, score) in enumerate(top30):
            print(f"  [{i+1}/{len(top30)}] {name}...", end=" ", flush=True)
            ok = run_waltz(seq)
            if ok:
                passed.append((name, seq, score))
                print("PASS")
            elif ok is False:
                print("FAIL")
            else:
                print("ERR")
        
        # Output
        out_fa = os.path.join(OUT_ROOT, f"{cand_id}_deA_WALTZ_passed.fa")
        out_all = os.path.join(OUT_ROOT, f"{cand_id}_deA_all.fa")
        with open(out_fa, "w") as f:
            for name, seq, score in passed:
                muts = sum(1 for a, b in zip(cfg["seq_orig"], seq) if a != b)
                a_new = seq.count("A")
                f.write(f">{name} muts={muts} A={a_new} agg_score={score}\n{seq}\n")
        with open(out_all, "w") as f:
            for name, seq, score in scored[:50]:
                f.write(f">{name} score={score}\n{seq}\n")
        print(f"\n  OUTPUT: {out_fa} ({len(passed)} WALTZ passed)")
        print(f"  ALL top50: {out_all}")
        
        # Summary
        if passed:
            best = passed[0]
            a_orig = cfg["seq_orig"].count("A")
            a_new = best[1].count("A")
            print(f"  BEST: {best[0]} A: {a_orig}→{a_new} ({(a_orig-a_new)/a_orig*100:.0f}% reduction)")

if __name__ == "__main__":
    main()
