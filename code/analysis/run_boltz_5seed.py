import os, subprocess, shutil

# Paths are configured through the environment so that the script is portable.
BASE = os.environ.get("DEA_REDESIGN_DIR", os.path.join("results", "deA_redesign_1_2"))
BOLTZ = os.environ.get("BOLTZ_BIN", "boltz")
IN_ROOT = os.path.join(BASE, "boltz_pdb_fullVEGF")
OUT_B2 = os.path.join(BASE, "boltz2_fullVEGF_test")

# new seeds to add (we already have 0,11,42 for both models)
NEW_SEEDS = (7, 23)
CANDS = ("1star", "2star", "4star")

log = os.path.join(OUT_B2, "run_boltz_5seed.log")
with open(log, "w", encoding="utf-8") as L:
    for cand in CANDS:
        # copy existing s0 yaml as the template for new seeds (same sequence, seed via CLI)
        src_yaml = os.path.join(IN_ROOT, f"{cand}_fullVEGF_s0", "inputs", f"{cand}_fullVEGF_s0.yaml")
        with open(src_yaml) as f:
            ytext = f.read()
        for seed in NEW_SEEDS:
            # --- Boltz-1 output into boltz_pdb_fullVEGF/{cand}_fullVEGF_s{seed} (mirror existing) ---
            b1_in = os.path.join(IN_ROOT, f"{cand}_fullVEGF_s{seed}", "inputs")
            os.makedirs(b1_in, exist_ok=True)
            b1_yaml = os.path.join(b1_in, f"{cand}_fullVEGF_s{seed}.yaml")
            with open(b1_yaml, "w") as f:
                f.write(ytext)
            b1_out = os.path.join(IN_ROOT, f"{cand}_fullVEGF_s{seed}", "out")
            os.makedirs(b1_out, exist_ok=True)
            # --- Boltz-2 output ---
            b2_out = os.path.join(OUT_B2, cand, f"boltz_results_{cand}_fullVEGF_s{seed}")
            os.makedirs(b2_out, exist_ok=True)

            for model, outd in (("boltz1", b1_out), ("boltz2", b2_out)):
                cmd = [
                    BOLTZ, "predict", "--model", model, "--no_kernels",
                    "--write_full_pde", "--diffusion_samples", "1", "--recycling_steps", "3",
                    "--seed", str(seed), "--out_dir", outd, b1_yaml,
                ]
                L.write(f"\n### {cand} s{seed} {model}\n")
                L.flush()
                r = subprocess.run(cmd, capture_output=True, text=True)
                L.write(f"EXIT={r.returncode}\n")
                L.write((r.stdout or "")[-800:] + "\n")
                L.write((r.stderr or "")[-800:] + "\n")
                L.flush()
                print(f"{cand} s{seed} {model} -> EXIT {r.returncode}")
    L.write("\nALL DONE\n")
print("5-SEED GRID COMPLETE")
