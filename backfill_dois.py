# -*- coding: utf-8 -*-
"""
backfill_dois.py — fill CSBJ availability placeholders after GitHub + Zenodo upload.

Usage (dry-run, prints what would change):
    python backfill_dois.py \
        --repo-url https://github.com/yangsitaoasprin/anti-vegf-ai \
        --github-doi 10.5281/zenodo.1234567 \
        --zenodo-doi 10.5281/zenodo.7654321

Apply the changes:
    python backfill_dois.py ... --apply

NOTE on the two DOIs:
  * --github-doi : the Zenodo-issued archive DOI of the *code* repo
                   (created via a GitHub "release" -> auto-minted). It replaces
                   `10.5281/zenodo.XXXXXXX` in the manuscripts.
  * --zenodo-doi : the DOI of the *raw 311 MB model outputs* deposit on Zenodo.
                   It replaces `10.5281/zenodo.YYYYYYYY` (manuscripts) and
                   `10.5281/zenodo.XXXXXXX` (results_manifest.md).
  The token strings differ by file, so replacements are file-scoped (see RULES).
"""

import argparse
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))  # anti_VEGF_AI_project/


def p(*rel):
    return os.path.join(ROOT, *rel)


# (file, [(old_substring, new_substring), ...])  -- file-scoped replacements
RULES = [
    # --- English manuscript: XXXXXXX=GitHub code DOI, YYYYYYYY=Zenodo raw DOI ---
    (p("论文", "Manuscript_v7.md"), [
        ("10.5281/zenodo.XXXXXXX", "{github_doi}"),
        ("10.5281/zenodo.YYYYYYYY", "{zenodo_doi}"),
    ]),
    # --- Chinese manuscript: same mapping ---
    (p("论文", "Manuscript_ZH.md"), [
        ("10.5281/zenodo.XXXXXXX", "{github_doi}"),
        ("10.5281/zenodo.YYYYYYYY", "{zenodo_doi}"),
    ]),
    # --- results_manifest.md: XXXXXXX = Zenodo raw DOI (NOT the github one) ---
    (p("repro", "results_manifest.md"), [
        ("10.5281/zenodo.XXXXXXX", "{zenodo_doi}"),
    ]),
    # --- CODE_AUTHORS.md: repo URL + zenodo raw DOI ---
    (p("repro", "CODE_AUTHORS.md"), [
        ("`<repo-URL>`", "`<{repo_url}>`"),
        ("`<zenodo-DOI>`", "`<{zenodo_doi}>`"),
    ]),
    # --- README.md: two descriptive placeholders -> concrete zenodo DOI ---
    (p("repro", "README.md"), [
        ("(DOI placeholder in `results_manifest.md`)",
         "(DOI: `https://doi.org/{zenodo_doi}`; see `results_manifest.md`)"),
        ("(DOI to be added after deposit)", "(`https://doi.org/{zenodo_doi}`)"),
    ]),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-url", required=True, help="public GitHub repo URL")
    ap.add_argument("--github-doi", required=True, help="code repo archive DOI (10.5281/zenodo.NNNN)")
    ap.add_argument("--zenodo-doi", required=True, help="raw-data Zenodo DOI (10.5281/zenodo.NNNN)")
    ap.add_argument("--apply", action="store_true", help="actually write files (default: dry-run)")
    args = ap.parse_args()

    fmt = {
        "repo_url": args.repo_url.rstrip("/"),
        "github_doi": args.github_doi.strip(),
        "zenodo_doi": args.zenodo_doi.strip(),
    }

    print("=== backfill DOIs ({} mode) ===".format("APPLY" if args.apply else "DRY-RUN"))
    print("  repo-url  :", fmt["repo_url"])
    print("  github-doi:", fmt["github_doi"])
    print("  zenodo-doi:", fmt["zenodo_doi"])
    print()

    changed = []
    for fpath, subs in RULES:
        if not os.path.isfile(fpath):
            print("[skip] not found:", fpath)
            continue
        with open(fpath, encoding="utf-8") as f:
            txt = f.read()
        original = txt
        for old, new in subs:
            new = new.format(**fmt)
            if old in txt:
                txt = txt.replace(old, new)
                print("[replace] {}  ::  {!r} -> {!r}".format(os.path.relpath(fpath, ROOT), old, new))
            else:
                print("[no-match] {}  ::  {!r}".format(os.path.relpath(fpath, ROOT), old))
        if txt != original:
            if args.apply:
                with open(fpath, "w", encoding="utf-8") as f:
                    f.write(txt)
                print("          wrote", os.path.relpath(fpath, ROOT))
            changed.append(os.path.relpath(fpath, ROOT))
        print()

    if changed:
        if args.apply:
            # Show what changed in the working tree
            print("=== git diff --stat (changes written) ===")
            try:
                out = subprocess.run(["git", "-C", ROOT, "diff", "--stat"],
                                     capture_output=True, text=True)
                print(out.stdout.strip() or "(not a git repo, or nothing staged)")
            except Exception:
                pass
            print("\nNext: commit these backfilled files, then push:")
            print("  git add " + " ".join(changed))
            print('  git commit -m "Backfill GitHub + Zenodo DOIs"')
        else:
            print("The above files WOULD change. Preview with: git diff")
            print("Apply with:  python backfill_dois.py ... --apply")
    else:
        print("Nothing to change (all placeholders already filled or absent).")

    print("Done. Re-run without --apply to preview, with --apply to write.")


if __name__ == "__main__":
    main()
