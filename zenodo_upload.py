#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zenodo_upload.py - stage and publish the two Zenodo data deposits for the manuscript.

WHY THIS EXISTS
    The agent process on this machine cannot resolve zenodo.org (DNS blocked), so the
    two uploads have to run from YOUR terminal.  This script turns ~15 minutes of
    drag-and-drop and metadata typing into one command, and it reads the resulting DOIs
    straight from the API - so no DOI is ever retyped by hand.

TWO-PHASE BY DESIGN
    Zenodo freezes metadata the moment a record is published (a fix costs a new version
    with a NEW DOI).  So:
        run 1 (no flags)          -> creates the two drafts, uploads the files, sets
                                     metadata, PRINTS the draft URLs, publishes nothing
        go look at them in the browser, fix anything you dislike
        run 2 (--publish)         -> publishes the staged drafts and prints the DOIs

USAGE
    # get a token: https://zenodo.org/account/settings/applications/tokens/new/
    #   scope: deposit:write   (tick "deposit:actions" too if offered)

    export ZENODO_TOKEN=xxxxxxxxxxxxxxxxxxxxx          # or pass --token
    python zenodo_upload.py                            # stage (safe, no publish)
    python zenodo_upload.py --publish                   # publish + print DOIs
    python zenodo_upload.py --publish --emit dois.txt   # also write them to a file

    # dry rehearsal against Zenodo's TEST instance (mints 10.5072/... DOIs, harmless):
    python zenodo_upload.py --sandbox --token <sandbox-token>

    # if a previous run left drafts behind and you want to reuse them instead of
    # creating duplicates:
    python zenodo_upload.py --use 1234567 7654321      # publish these ids as-is

NOTES
    * --sandbox uses https://sandbox.zenodo.org ; DOIs minted there do NOT resolve
      publicly and must never go into the manuscript.
    * The script prints only a masked token prefix, never the token itself.
    * Uploads stream the file with an explicit Content-Length (no full-file buffering).
"""
import argparse
import datetime
import http.client
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
DEPO = os.path.join(ROOT, "zenodo_deposit")

REPO_URL = "https://github.com/yangsitaoasprin/anti-vegf-ai"
MANUSCRIPT_TITLE = ("Developability-gated virtual screening of de novo "
                    "anti-VEGF-A miniprotein binders")
CREATOR = {
    "name": "Yang, Sitao",
    "affiliation": "The Third Affiliated Hospital of Dali University, Dali, Yunnan, China; "
                   "School of Pharmaceutical Science, Dali University, Dali, Yunnan, China",
    "orcid": "0009-0006-6442-4250",
}
LICENSE = "cc-by-4.0"

DEPOSITS = [
    {
        "key": "raw",
        "filename": "anti_VEGF_AI_raw_predictions_v1.zip",
        "title": 'Raw model outputs for "%s"' % MANUSCRIPT_TITLE,
        "keywords": ["VEGF-A", "miniprotein binder", "Boltz-1", "OpenFold3",
                     "MM-GBSA", "de novo protein design"],
        "description": (
            "<p>Complete raw structure-prediction outputs underlying the manuscript "
            "<em>&ldquo;%s&rdquo;</em>.</p>"
            "<p>This archive holds the full <b>Boltz-1</b> predicted complexes and the "
            "complete <b>OpenFold3</b> predictions (mmCIF structures with their "
            "per-prediction confidence JSON) for every design candidate and every seed, "
            "together with the Schr&ouml;dinger <b>MM-GBSA</b> working directories used "
            "for the per-pose energy decomposition. Every numerical value reported in the "
            "manuscript is derived from these files.</p>"
            "<p>A manifest of the archive, the final candidate table, the per-seed "
            "MM-GBSA decomposition and all figure-generation code are in the companion "
            "GitHub repository. Code is released under MIT; these data under CC-BY-4.0.</p>"
            % MANUSCRIPT_TITLE
        ),
    },
    {
        "key": "apr",
        "filename": "anti_VEGF_AI_apr_ccs_recovery_v1.zip",
        "title": ('APR-Score and CCS recovery - analysis record for the manuscript '
                  '"%s"' % MANUSCRIPT_TITLE),
        "keywords": ["amyloid", "aggregation propensity", "APR-Score", "CCS",
                     "protein design", "reproducibility"],
        "description": (
            "<p>Analysis record behind the APR-Score / CCS aggregation scoring and "
            "Supplementary Fig. S1 of the manuscript <em>&ldquo;%s&rdquo;</em>.</p>"
            "<p>The archive contains the recovery scripts run against the published "
            "APR-Score binary, the derived tables they regenerate, the verification "
            "ladder used to test the reconstruction, and the internal reports that "
            "document both how the official executable fails as distributed and how the "
            "17-of-18-dimensional featuriser was reconstructed from published scale "
            "descriptions. Interpretability of the one non-identifiable dimension, and "
            "its effect on the ranking reported in the paper, are stated explicitly in "
            "the accompanying report.</p>"
            "<p>Code is released under MIT; these data under CC-BY-4.0.</p>"
            % MANUSCRIPT_TITLE
        ),
    },
]


def mask(tok):
    return (tok[:4] + "***" + tok[-4:]) if tok and len(tok) > 12 else "(short)"


class Zenodo(object):
    def __init__(self, token, sandbox=False, verbose=True):
        self.base = "https://sandbox.zenodo.org" if sandbox else "https://zenodo.org"
        self.token = token
        self.verbose = verbose

    def _req(self, path, method="GET", payload=None, timeout=120):
        url = self.base + path
        data = json.dumps(payload).encode() if payload is not None else None
        req = urllib.request.Request(url, data=data, method=method, headers={
            "Authorization": "Bearer " + self.token,
            "Content-Type": "application/json",
            "User-Agent": "anti-vegf-ai-deposit/1.0",
        })
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                body = r.read().decode()
                return r.status, (json.loads(body) if body else {})
        except urllib.error.HTTPError as e:
            raw = e.read().decode(errors="replace")
            try:
                raw = json.loads(raw)
            except Exception:
                pass
            return e.code, raw

    def check(self):
        return self._req("/api/depositions?size=1")

    def create(self):
        return self._req("/api/depositions", "POST", {})

    def set_metadata(self, dep_id, meta):
        return self._req("/api/depositions/%s" % dep_id, "PUT", {"metadata": meta})

    def publish(self, dep_id):
        return self._req("/api/depositions/%s/actions/publish" % dep_id, "POST", timeout=300)

    def upload(self, bucket, local_path):
        name = os.path.basename(local_path)
        size = os.path.getsize(local_path)
        u = urllib.parse.urlsplit(bucket.rstrip("/") + "/" + urllib.parse.quote(name))
        conn = http.client.HTTPSConnection(u.hostname, u.port or 443, timeout=1800)
        path = u.path
        fh = open(local_path, "rb")
        try:
            conn.putrequest("PUT", path, skip_accept_encoding=True)
            conn.putheader("Authorization", "Bearer " + self.token)
            conn.putheader("Content-Type", "application/octet-stream")
            conn.putheader("Content-Length", str(size))
            conn.endheaders()
            sent = 0
            mark = 0
            while True:
                chunk = fh.read(1 << 20)
                if not chunk:
                    break
                conn.send(chunk)
                sent += len(chunk)
                if self.verbose and size >= (5 << 20) and sent - mark >= (10 << 20):
                    mark = sent
                    sys.stderr.write("      ... %d / %d MB (%.0f%%)\n"
                                     % (sent >> 20, size >> 20, 100.0 * sent / size))
            resp = conn.getresponse()
            body = resp.read().decode(errors="replace")
            return resp.status, body
        finally:
            fh.close()
            conn.close()


def meta_for(d, today):
    return {
        "title": d["title"],
        "upload_type": "dataset",
        "publication_date": today,
        "description": d["description"],
        "creators": [CREATOR],
        "access_right": "open",
        "license": LICENSE,
        "keywords": d["keywords"],
        "related_identifiers": [
            {"identifier": REPO_URL, "relation": "isSupplementedBy",
             "resource_type": "software"},
        ],
        "notes": ("This record is one of three deposits supporting the manuscript: the "
                  "code repository (GitHub, archived by Zenodo), the raw model outputs, "
                  "and this analysis record."),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--token", default=os.environ.get("ZENODO_TOKEN", ""),
                    help="Zenodo personal access token (or set ZENODO_TOKEN)")
    ap.add_argument("--sandbox", action="store_true", help="use sandbox.zenodo.org (test DOIs)")
    ap.add_argument("--publish", action="store_true", help="publish the staged drafts")
    ap.add_argument("--use", nargs="*", default=[], metavar="ID",
                    help="operate on these existing deposition ids instead of creating new")
    ap.add_argument("--emit", default="", metavar="FILE", help="write the DOIs to FILE")
    ap.add_argument("--force-large", action="store_true",
                    help="proceed even if a zip is above the 50 GB per-file limit")
    args = ap.parse_args()

    if not args.token:
        sys.exit("ERROR: no token. Create one at "
                 f"{'https://sandbox.zenodo.org' if args.sandbox else 'https://zenodo.org'}"
                 "/account/settings/applications/tokens/new/ (scope: deposit:write), then\n"
                 "       export ZENODO_TOKEN=...   or pass --token ...")

    print("=== Zenodo deposit: %s ===" % ("SANDBOX (test DOIs - do NOT use in the paper)"
                                          if args.sandbox else "PRODUCTION"))
    print("  token      : %s" % mask(args.token))
    print("  files dir  : %s" % DEPO)
    print("  mode       : %s" % ("PUBLISH" if args.publish else "STAGE ONLY (nothing published)"))
    print()

    z = Zenodo(args.token, sandbox=args.sandbox)
    st, _ = z.check()
    if st != 200:
        sys.exit("ERROR: token rejected (%d). Check scope deposit:write." % st)
    print("  token accepted.\n")

    today = datetime.date.today().isoformat()
    results = []

    # ---------------------------------------------------------------- existing ids?
    if args.use:
        todo = []
        for i, dep_id in enumerate(args.use):
            st, dep = z._req("/api/depositions/%s" % dep_id)
            if st != 200:
                sys.exit("ERROR: cannot read deposition %s (%d)" % (dep_id, st))
            todo.append((DEPOSITS[i] if i < len(DEPOSITS) else DEPOSITS[-1], dep))
    else:
        todo = []
        for d in DEPOSITS:
            path = os.path.join(DEPO, d["filename"])
            if not os.path.isfile(path):
                sys.exit("ERROR: missing file %s" % path)
            size = os.path.getsize(path)
            if size > 50 * (1 << 30) and not args.force_large:
                sys.exit("ERROR: %s is %.1f GB, above Zenodo's 50 GB per-file limit"
                         % (d["filename"], size / float(1 << 30)))
            print("--- %s (%s, %.1f MB) ---" % (d["key"], d["filename"], size / 1048576.0))
            st, dep = z.create()
            if st not in (200, 201):
                sys.exit("ERROR: create failed (%d): %s" % (st, json.dumps(dep)[:400]))
            dep_id = dep["id"]
            bucket = dep["links"]["bucket"]
            pre = (dep.get("metadata") or {}).get("prereserve_doi") or {}
            print("    created draft id %s" % dep_id)
            if pre:
                print("    reserved DOI  : %s" % pre.get("doi"))
            print("    uploading ...")
            st, body = z.upload(bucket, path)
            if st not in (200, 201):
                sys.exit("ERROR: upload of %s failed (%d): %s"
                         % (d["filename"], st, str(body)[:400]))
            print("    uploaded ok")
            st, res = z.set_metadata(dep_id, meta_for(d, today))
            if st not in (200, 201):
                sys.exit("ERROR: metadata rejected (%d): %s" % (st, json.dumps(res)[:600]))
            print("    metadata set")
            todo.append((d, res))
            print()

    # ---------------------------------------------------------------- publish / report
    print("=== result ===")
    for d, dep in todo:
        dep_id = dep["id"]
        html = (dep.get("links") or {}).get("html", "")
        if args.publish:
            st, res = z.publish(dep_id)
            if st not in (200, 201, 202):
                print("  [FAIL] %-4s publish rejected (%d): %s"
                      % (d["key"], st, json.dumps(res)[:400]))
                continue
            doi = res.get("doi") or (res.get("metadata") or {}).get("doi", "")
            print("  [LIVE] %-4s %s" % (d["key"], doi))
            print("         %s" % ((res.get("links") or {}).get("html", html)))
            results.append((d["key"], doi))
        else:
            pre = (dep.get("metadata") or {}).get("prereserve_doi") or {}
            print("  [DRAFT] %-4s id=%s" % (d["key"], dep_id))
            print("          %s" % html)
            if pre:
                print("          reserved DOI: %s" % pre.get("doi"))
    print()

    if not args.publish:
        print("Nothing was published. Review the drafts in the browser, then re-run:")
        print("  python zenodo_upload.py --publish --use %s"
              % " ".join(str(dep["id"]) for _, dep in todo))
        print("Metadata is FROZEN at publish time - check the titles against the manuscript.")
    else:
        if args.emit and results:
            with open(args.emit, "w", encoding="utf-8") as f:
                for k, doi in results:
                    f.write("%s\t%s\n" % (k, doi))
            print("DOIs written to %s" % args.emit)
        print("Feed these to backfill_dois.py as --zenodo-doi (raw) and --apr-doi (record).")


if __name__ == "__main__":
    main()
