#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zenodo_deposit.py — 用 Zenodo 新版 InvenioRDM REST API 完成三件事：

  1) --stage          建草稿 + 流式上传 zip + 写元数据（**不发布**）
  2) --publish       把已 stage 的草稿发布，拿版本 DOI
  3) --edit-record   订正一条**已发布**记录的元数据（DOI 不变、不产生新版本）
  另 --discard 删草稿；--list 列出我的 records。

为什么不用旧脚本 `zenodo_upload.py`：它打的是 legacy `/api/depositions`，
该路径在生产上已 **404**（真实前缀是 `/api/deposit/depositions`）；
新版 `POST /api/records` 才是官方支持路径。载荷字段名也整代不同（见下 METADATA 映射）。

DNS 兜底：本机对 `zenodo.org` / `api.zenodo.org` 的解析被投毒（Errno 11004/11001），
故内置 PinnedHTTPSConnection —— 连已知 IP，但 SNI 与 Host 仍写 `zenodo.org`，
证书 `*.zenodo.org` 照常校验通过。默认自动尝试：先普通连接，失败即切 pin。

用法
----
  export ZENODO_TOKEN=...
  python zenodo_deposit.py --stage                 # 看它建了什么草稿（不发布）
  python zenodo_deposit.py --publish --use 123 456 # 发布
  python zenodo_deposit.py --edit-record 22866149  # 订正已发布记录元数据
"""
import argparse
import datetime
import http.client
import json
import os
import socket
import ssl
import sys
import urllib.error
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
DEPO = os.path.join(ROOT, "zenodo_deposit")

PROD_HOST = "zenodo.org"
SANDBOX_HOST = "sandbox.zenodo.org"
# 实测可用（2026-09-21，SNI=zenodo.org 均通过 *.zenodo.org 证书校验）
PIN_IPS = ["137.138.153.219", "137.138.52.235", "188.184.103.118", "188.184.98.114"]

REPO_URL = "https://github.com/yangsitaoasprin/anti-vegf-ai"
CODE_DOI = "10.5281/zenodo.22866148"          # 代码仓概念 DOI
CODE_DOI_V12 = "10.5281/zenodo.22866393"      # 代码仓 v1.2 版本 DOI
MANUSCRIPT_TITLE = ("Developability-gated virtual screening of de novo "
                    "anti-VEGF-A miniprotein binders")
CREATOR = {
    "person_or_org": {
        "type": "personal",
        "family_name": "Yang",
        "given_name": "Sitao",
        "identifiers": [{"identifier": "0009-0006-6442-4250", "scheme": "orcid"}],
    },
    "affiliations": [{"name": "The Third Affiliated Hospital of Dali University, "
                              "Dali, Yunnan, China; College of Pharmacy, "
                              "Dali University, Dali, Yunnan, China"}],
}

DEPOSITS = [
    {
        "key": "raw",
        "filename": "anti_VEGF_AI_raw_predictions_v1.zip",
        "title": 'Raw model outputs for "%s"' % MANUSCRIPT_TITLE,
        "subjects": ["VEGF-A", "miniprotein binder", "Boltz-1", "OpenFold3",
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
            "<p>This record is one of three deposits supporting the manuscript: the code "
            "repository (GitHub, archived by Zenodo), the raw model outputs, and the "
            "APR-Score / CCS analysis record.</p>"
            % MANUSCRIPT_TITLE
        ),
    },
    {
        "key": "apr",
        "filename": "anti_VEGF_AI_apr_ccs_recovery_v1.zip",
        "title": ('APR-Score and CCS recovery - analysis record for the manuscript '
                  '"%s"' % MANUSCRIPT_TITLE),
        "subjects": ["amyloid", "aggregation propensity", "APR-Score", "CCS",
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
            "<p>This record is one of three deposits supporting the manuscript: the code "
            "repository (GitHub, archived by Zenodo), the raw model outputs, and the "
            "APR-Score / CCS analysis record.</p>"
            % MANUSCRIPT_TITLE
        ),
    },
]


def meta_new(d, today):
    """新版 InvenioRDM 元数据载荷。"""
    return {
        "title": d["title"],
        "publication_date": today,
        "resource_type": {"id": "dataset"},
        "creators": [CREATOR],
        "description": d["description"],
        "rights": [{"id": "cc-by-4.0"}],
        "subjects": [{"subject": s} for s in d["subjects"]],
        "version": "v1.0",
        "related_identifiers": [
            {"identifier": CODE_DOI, "scheme": "doi",
             "relation_type": {"id": "issupplementedby"},
             "resource_type": {"id": "software"}},
            {"identifier": REPO_URL, "scheme": "url",
             "relation_type": {"id": "issupplementedby"},
             "resource_type": {"id": "software"}},
        ],
    }


# ---------------------------------------------------------------- HTTP with DNS fallback
class PinnedHTTPSConnection(http.client.HTTPSConnection):
    def __init__(self, host, ip, port=443, timeout=60, **kw):
        self._pin_ip = ip
        super().__init__(host, port, timeout=timeout, **kw)

    def connect(self):
        sock = socket.create_connection((self._pin_ip, self.port), self.timeout)
        ctx = self._context or ssl.create_default_context()
        self.sock = ctx.wrap_socket(sock, server_hostname=self.host)


class ProxyHTTPSConnection(http.client.HTTPSConnection):
    """经本地 HTTP 代理的 CONNECT 隧道连目标，SNI/Host 仍是真域名。

    为什么需要它（2026-09-21 实测）：裸连（pin）到 CERN 只有 ~20 KB/s，
    102 MB 的包传 10-40% 就被重置，四个 IP 全失败。而本机系统代理
    127.0.0.1:5780 能 CONNECT zenodo.org 并给到 120-203 KB/s，
    约 9 分钟就能传完。出口依然偶发被重置，故调用方要多试几轮。
    """

    def __init__(self, host, proxy_host, proxy_port, port=443, timeout=60, **kw):
        self._ph = proxy_host
        self._pp = int(proxy_port)
        super().__init__(host, port, timeout=timeout, **kw)

    def connect(self):
        sock = socket.create_connection((self._ph, self._pp), self.timeout)
        try:
            sock.sendall(("CONNECT %s:%d HTTP/1.1\r\nHost: %s:%d\r\n\r\n"
                          % (self.host, self.port, self.host, self.port)).encode())
            resp = b""
            while b"\r\n\r\n" not in resp:
                b = sock.recv(4096)
                if not b:
                    raise OSError("proxy closed the connection during CONNECT")
                resp += b
            first = resp.split(b"\r\n")[0].decode(errors="replace")
            if " 200" not in first:
                raise OSError("proxy refused CONNECT: %s" % first)
        except Exception:
            sock.close()
            raise
        ctx = self._context or ssl.create_default_context()
        self.sock = ctx.wrap_socket(sock, server_hostname=self.host)


class Client(object):
    """所有请求都经这里。mode: auto | direct | pin；proxy=(host, port) 时优先走代理。"""

    def __init__(self, host, token, mode="auto", verbose=True, proxy=None,
                 proxy_tries=8):
        self.host = host
        self.token = token
        self.mode = mode
        self.verbose = verbose
        self.chosen = None
        self.proxy = proxy
        self.proxy_tries = proxy_tries

    # ---- 低层：经本地代理隧道
    def _via_proxy(self, method, path, payload=None, timeout=120, tries=3):
        last = None
        for _ in range(tries):
            try:
                c = ProxyHTTPSConnection(self.host, self.proxy[0], self.proxy[1],
                                         timeout=timeout)
                body = json.dumps(payload).encode() if payload is not None else None
                c.request(method, path, body=body, headers={
                    "Authorization": "Bearer " + self.token,
                    "Content-Type": "application/json",
                    "User-Agent": "anti-vegf-ai-deposit/2.0"})
                r = c.getresponse()
                raw = r.read().decode(errors="replace")
                st = r.status
                c.close()
                try:
                    return st, (json.loads(raw) if raw else {})
                except Exception:
                    return st, raw
            except Exception as e:
                last = e
        raise last

    # ---- 低层：普通 urllib（direct）
    def _direct(self, method, path, payload=None, timeout=120):
        url = "https://%s%s" % (self.host, path)
        data = json.dumps(payload).encode() if payload is not None else None
        req = urllib.request.Request(url, data=data, method=method, headers={
            "Authorization": "Bearer " + self.token,
            "Content-Type": "application/json",
            "User-Agent": "anti-vegf-ai-deposit/2.0",
        })
        op = urllib.request.build_opener(urllib.request.ProxyHandler({}))
        try:
            with op.open(req, timeout=timeout) as r:
                body = r.read().decode()
                return r.status, (json.loads(body) if body else {})
        except urllib.error.HTTPError as e:
            raw = e.read().decode(errors="replace")
            try:
                raw = json.loads(raw)
            except Exception:
                pass
            return e.code, raw

    # ---- 低层：pin 连接（自建 socket，SNI=真实域名）
    def _pin(self, method, path, payload=None, timeout=120):
        last = None
        ips = PIN_IPS if self.host == PROD_HOST else [None]
        for ip in ips:
            if ip is None:
                return self._direct(method, path, payload, timeout)
            try:
                c = PinnedHTTPSConnection(self.host, ip, timeout=timeout)
                body = json.dumps(payload).encode() if payload is not None else None
                hdrs = {"Authorization": "Bearer " + self.token,
                        "Content-Type": "application/json",
                        "User-Agent": "anti-vegf-ai-deposit/2.0"}
                c.request(method, path, body=body, headers=hdrs)
                r = c.getresponse()
                raw = r.read().decode(errors="replace")
                st = r.status
                c.close()
                if self.verbose and self.chosen != ip:
                    self.chosen = ip
                try:
                    return st, (json.loads(raw) if raw else {})
                except Exception:
                    return st, raw
            except Exception as e:
                last = e
        raise last

    def req(self, method, path, payload=None, timeout=120):
        if self.proxy:
            return self._via_proxy(method, path, payload, timeout)
        if self.mode == "direct":
            return self._direct(method, path, payload, timeout)
        if self.mode == "pin":
            return self._pin(method, path, payload, timeout)
        # auto: 先 direct，DNS 挂了再 pin
        try:
            return self._direct(method, path, payload, timeout)
        except Exception as e:
            if self.verbose:
                sys.stderr.write("      (direct 失败 %s，切 IP 直连)\n" % type(e).__name__)
            self.mode = "pin"
            return self._pin(method, path, payload, timeout)

    # ---- 大文件流式 PUT（分块，显式 Content-Length），带 pin 兜底
    def put_file(self, path, local_path):
        name = os.path.basename(local_path)
        size = os.path.getsize(local_path)
        if self.proxy:
            # 经本地代理：实测比 pin 快 6-10 倍，但出口仍会偶发重置。
            # API 不支持断点续传，所以每轮都从头开始 —— 多给几轮机会。
            routes = [("proxy", None)] * max(1, self.proxy_tries)
        elif self.mode == "direct":
            routes = [("direct", None)]
        else:
            routes = [("pin", ip) for ip in
                      (PIN_IPS if self.host == PROD_HOST else [None])]
        last = None
        for attempt, route in enumerate(routes, 1):
            fh = open(local_path, "rb")
            conn = None
            try:
                if route[0] == "proxy":
                    conn = ProxyHTTPSConnection(self.host, self.proxy[0], self.proxy[1],
                                                timeout=1800)
                elif route[0] == "direct":
                    conn = http.client.HTTPSConnection(self.host, 443, timeout=1800)
                else:
                    conn = PinnedHTTPSConnection(self.host, route[1], timeout=1800)
                if attempt > 1 and self.verbose:
                    sys.stderr.write("      (第 %d 轮重试，从头开始：%s)\n"
                                     % (attempt, route[0]))
                # 注意：putrequest 默认已带 Host（取 self.host），不要再手动加一个，
                # 否则会发出两个 Host 头。
                # ⚠️ 必须带 User-Agent：实测 Zenodo 的 WAF 对「没有 UA」和「伪装浏览器的
                # UA（Mozilla/...）」都回 403 HTML 拦截页 —— 只有脚本型 UA 放行。
                conn.putrequest("PUT", path, skip_accept_encoding=True)
                conn.putheader("User-Agent", "anti-vegf-ai-deposit/2.0")
                conn.putheader("Authorization", "Bearer " + self.token)
                conn.putheader("Content-Type", "application/octet-stream")
                conn.putheader("Content-Length", str(size))
                conn.endheaders()
                sent, mark = 0, 0
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
                try:
                    resp = conn.getresponse()
                    raw = resp.read().decode(errors="replace")
                    return resp.status, raw
                except Exception as exc:
                    return 0, ("connection aborted before the server answered (%s: %s). "
                               "Most likely the token lacks deposit:write, or the network "
                               "dropped. The draft exists and is empty - re-run with "
                               "--use <draft id>." % (type(exc).__name__, exc))
            except Exception as e:
                last = e
                sys.stderr.write("      (上传走 %s 失败：%s)\n" % (route[0], trunc(e)))
            finally:
                fh.close()
                if conn is not None:
                    try:
                        conn.close()
                    except Exception:
                        pass
        raise last if last else RuntimeError("upload failed")


def trunc(e, n=140):
    s = "%s: %s" % (type(e).__name__, e)
    return s if len(s) <= n else s[:n] + "…"


def mask(t):
    return (t[:4] + "***" + t[-4:]) if t and len(t) > 12 else "(short)"


# ---------------------------------------------------------------- operations
def stage(z, today, only=None):
    out = []
    for d in DEPOSITS:
        if only and d["key"] not in only:
            continue
        path = os.path.join(DEPO, d["filename"])
        if not os.path.isfile(path):
            sys.exit("ERROR: 缺文件 %s" % path)
        size = os.path.getsize(path)
        print("--- %s (%s, %.1f MB) ---" % (d["key"], d["filename"], size / 1048576.0))
        payload = {"access": {"record": "public", "files": "public"},
                   "files": {"enabled": True},
                   "metadata": meta_new(d, today)}
        st, dep = z.req("POST", "/api/records", payload)
        if st not in (200, 201):
            sys.exit("ERROR: 建草稿失败 (%s): %s" % (st, json.dumps(dep)[:500]))
        rid = dep["id"]
        links = dep.get("links") or {}
        print("    草稿 id      : %s" % rid)
        print("    草稿页面     : %s" % links.get("self_html", ""))
        pids = dep.get("pids") or {}
        if pids.get("doi"):
            print("    预留 DOI     : %s" % pids["doi"].get("identifier"))
        # 1) init file
        st, r = z.req("POST", "/api/records/%s/draft/files" % rid,
                      [{"key": d["filename"]}])
        if st not in (200, 201):
            sys.exit("ERROR: 初始化文件失败 (%s): %s" % (st, json.dumps(r)[:400]))
        # 2) upload content
        up = links.get("files") or "/api/records/%s/draft/files" % rid
        up_path = urllib.parse.urlsplit(up).path
        target = "%s/%s/content" % (up_path.rstrip("/"),
                                    urllib.parse.quote(d["filename"]))
        print("    上传中 ...")
        st, body = z.put_file(target, path)
        if st not in (200, 201):
            sys.exit("ERROR: 上传失败 (%s): %s" % (st, str(body)[:400]))
        # 3) commit
        st, r = z.req("POST", "%s/%s/commit" % (up_path.rstrip("/"),
                                               urllib.parse.quote(d["filename"])))
        if st not in (200, 201):
            sys.exit("ERROR: commit 失败 (%s): %s" % (st, json.dumps(r)[:400]))
        print("    上传完成，服务端登记：size=%s checksum=%s"
              % ((r or {}).get("size"), (r or {}).get("checksum")))
        out.append((d, rid, (r or {}).get("size"), (r or {}).get("checksum")))
        print()
    return out


def _md5(path, chunk=1 << 20):
    import hashlib
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _norm_files(dep):
    """文件列表：兼容 RDM 的 {entries:{...}} 与 legacy 的 [...]。"""
    raw = dep.get("files")
    if isinstance(raw, dict):
        ent = raw.get("entries")
        if isinstance(ent, dict):
            return list(ent.values())
        if isinstance(ent, list):
            return ent
        return []
    return raw if isinstance(raw, list) else []


def _norm_meta(md):
    """把两种序列化（RDM / Zenodo legacy）归一，便于核对与打印。

    ⚠️ 实测：Zenodo 的 GET 返回的是 **legacy 形状**（creators[].name、
    license、keywords、related_identifiers[].relation），而写入要 **RDM 形状**
    （creators[].person_or_org、rights、subjects、relation_type）。两套必须都能读。
    """
    c0 = (md.get("creators") or [{}])[0]
    po = c0.get("person_or_org") or {}
    if po:
        name = ", ".join(x for x in (po.get("family_name"), po.get("given_name")) if x)
        affil = [a.get("name") for a in (c0.get("affiliations") or [])]
        orcid = [i.get("identifier") for i in (po.get("identifiers") or [])]
    else:
        name = c0.get("name")
        affil = [c0.get("affiliation")] if c0.get("affiliation") else []
        orcid = [c0.get("orcid")] if c0.get("orcid") else []
    lic = md.get("rights") or md.get("license") or []
    if isinstance(lic, dict):
        lic = [lic]
    lic = [x.get("id") if isinstance(x, dict) else x for x in lic]
    subj = md.get("subjects")
    if subj is None:
        subj = md.get("keywords") or []
    subj = [x.get("subject") if isinstance(x, dict) else x for x in subj]
    rel = [(x.get("relation_type", {}).get("id") or x.get("relation"), x.get("identifier"))
           for x in (md.get("related_identifiers") or [])]
    rt = md.get("resource_type") or {}
    return dict(title=md.get("title"), version=md.get("version"),
                publication_date=md.get("publication_date"),
                name=name, affil=affil, orcid=orcid, lic=lic, subj=subj, rel=rel,
                resource_type=(rt.get("id") or rt.get("type")))


def _show(tag, n):
    print("  [%s] id=%s" % (tag, n.get("_rid", "?")))
    print("      title    : %s" % (n["title"] or "")[:96])
    print("      creator  : %s   ORCID=%s" % (n["name"], n["orcid"]))
    print("      affil    : %s" % n["affil"])
    print("      license  : %s   type=%s" % (n["lic"], n["resource_type"]))
    print("      keywords : %d 个  %s" % (len(n["subj"]), n["subj"][:3]))
    print("      version  : %s   pubdate=%s" % (n["version"], n["publication_date"]))
    print("      related  : %s" % n["rel"])


def _sum_of(f):
    """checksum 可能是 'md5:xxx' 字符串，也可能是 {"algorithm":..,"checksum":..}。"""
    c = f.get("checksum")
    if isinstance(c, dict):
        return c.get("checksum")
    return (c or "").replace("md5:", "")


def verify(z, staged):
    print("=" * 92)
    print("草稿核对（发布前必看）—— 大小 + md5 双重比对")
    print("=" * 92)
    ok = True
    for d, rid, esize, esum in staged:
        st, dep = z.req("GET", "/api/records/%s/draft" % rid)
        if st not in (200, 201):
            print("  [%s] 读草稿失败 %s" % ((d["key"] if d else rid), st))
            ok = False
            continue
        md = dep.get("metadata") or {}
        n = _norm_meta(md)
        n["_rid"] = rid
        fl = _norm_files(dep)
        keys = [(f.get("key") or f.get("filename")) for f in fl]
        d2 = next((x for x in DEPOSITS if x["filename"] in keys), None)
        if d2 is None:
            print("  [%s] 草稿里的文件名认不出属于哪个沉积：%s" % (rid, keys))
            ok = False
            continue
        if d is not None and d2 is not d:
            print("      （按文件名识别为 %s，已纠正位置对应）" % d2["key"])
        d = d2
        # 期望值一律以本地磁盘为准。早先 --verify 按位置把 args.verify[i]
        # 配到 DEPOSITS[i]，只核对一个草稿时会把 apr 的 id 配上 raw 的
        # 102 MB 期望值，于是大小对不上、报假 FAIL。
        local = os.path.join(DEPO, d["filename"])
        if not os.path.isfile(local):
            print("  [%s] 本地文件不存在，无法比对：%s" % (d["key"], local))
            ok = False
            continue
        esize = os.path.getsize(local)
        emd5 = _md5(local)
        _show(d["key"], n)
        for f in fl:
            key = f.get("key") or f.get("filename")
            sz = f.get("size") if f.get("size") is not None else f.get("filesize")
            ck = _sum_of(f)
            s_ok = (sz == esize)
            m_ok = (ck == emd5)
            good = s_ok and m_ok
            ok = ok and good
            print("      file     : %s" % key)
            print("                 服务端 %s B  md5:%s" % (sz, ck))
            print("                 本地   %s B  md5:%s   -> %s"
                  % (esize, emd5, "PASS 一致" if good else "FAIL 不一致"))
        print()
    print("核对判定: %s" % ("PASS" if ok else "有问题，先别发布"))
    return ok


def publish(z, ids):
    for rid in ids:
        st, r = z.req("POST", "/api/records/%s/draft/actions/publish" % rid)
        if st not in (200, 201, 202):
            print("  [FAIL] %s 发布失败 (%s): %s" % (rid, st, json.dumps(r)[:400]))
            continue
        pids = r.get("pids") or {}
        doi = (pids.get("doi") or {}).get("identifier") or r.get("doi", "")
        print("  [LIVE] id=%s  DOI=%s" % (rid, doi))
        print("         %s" % ((r.get("links") or {}).get("self_html", "")))
    print()
    print("把 DOI 交给 backfill_dois.py：raw 传 --zenodo-doi，APR 传 --apr-doi。")


def discard(z, ids):
    for rid in ids:
        st, r = z.req("DELETE", "/api/records/%s/draft" % rid)
        print("  删除草稿 %s -> HTTP %s" % (rid, st))


KEYWORDS_CODE = ["VEGF-A", "miniprotein binder", "de novo protein design", "Boltz-1",
                 "OpenFold3", "MM-GBSA", "developability screening", "reproducibility"]


def legacy_to_new(md):
    """legacy 形状的 metadata → RDM 形状（写入 PUT/POST 用）。"""
    rt = md.get("resource_type") or {}
    rt_id = rt.get("id") or rt.get("type") or "software"
    lic = md.get("rights") or md.get("license") or []
    if isinstance(lic, dict):
        lic = [lic]
    lic_ids = [x.get("id") if isinstance(x, dict) else x for x in lic]
    subj = md.get("subjects")
    if subj is None:
        subj = md.get("keywords") or []
    subj = [x.get("subject") if isinstance(x, dict) else x for x in subj]
    c0 = (md.get("creators") or [{}])[0]
    if c0.get("person_or_org"):
        creators = [c0]
    else:
        nm = (c0.get("name") or "").strip()
        if "," in nm:
            fam, _, giv = nm.partition(",")
        else:
            parts = nm.split()
            fam, giv = (parts[-1] if parts else nm), " ".join(parts[:-1])
        ident = ([{"identifier": c0["orcid"], "scheme": "orcid"}]
                 if c0.get("orcid") else [])
        creators = [{"person_or_org": {"type": "personal",
                                       "family_name": fam.strip(),
                                       "given_name": giv.strip(),
                                       "identifiers": ident},
                     "affiliations": ([{"name": c0["affiliation"]}]
                                      if c0.get("affiliation") else [])}]
    rel = []
    for x in (md.get("related_identifiers") or []):
        ident = x.get("identifier") or ""
        r = ((x.get("relation_type") or {}).get("id") or x.get("relation") or "").lower()
        e = {"identifier": ident,
             "scheme": x.get("scheme") or ("doi" if ident.startswith("10.") else "url"),
             "relation_type": {"id": r}}
        rtv = x.get("resource_type")
        if isinstance(rtv, str):
            e["resource_type"] = {"id": rtv}
        elif isinstance(rtv, dict):
            e["resource_type"] = rtv
        rel.append(e)
    out = {"title": md.get("title"),
           "publication_date": md.get("publication_date"),
           "resource_type": {"id": rt_id},
           "creators": creators,
           "description": md.get("description"),
           "rights": [{"id": i} for i in lic_ids if i],
           "subjects": [{"subject": s} for s in subj if s],
           "related_identifiers": rel}
    if md.get("version"):
        out["version"] = md["version"]
    return out


def fix_record(z, rid, dry=True):
    """订正已发布记录的元数据：建 edit 草稿 -> 改 -> 发布。DOI 不变、不产生新版本。"""
    st, rec = z.req("GET", "/api/records/%s" % rid)
    if st != 200:
        sys.exit("ERROR: 读记录 %s 失败 (%s): %s" % (rid, st, json.dumps(rec)[:300]))
    before = rec.get("metadata") or {}
    nb = _norm_meta(before)
    nb["_rid"] = rid
    nb["title"] = "%s   (doi=%s)" % (nb["title"], rec.get("doi"))
    print("=" * 92)
    print("记录 %s 当前元数据" % rid)
    print("=" * 92)
    _show("before", nb)

    after = legacy_to_new(before)
    after["title"] = ('Reproducible code and supporting data for "%s"' % MANUSCRIPT_TITLE)
    after["creators"] = [CREATOR]
    after["rights"] = [{"id": "mit"}]
    after["subjects"] = [{"subject": s} for s in KEYWORDS_CODE]
    desc = after.get("description") or ""
    if "Licence: CC-BY-4.0." in desc:
        desc = desc.replace("Licence: CC-BY-4.0.", "Licence: MIT for the code.")
    after["description"] = desc

    na = _norm_meta(after)
    na["_rid"] = rid
    print()
    print("=" * 92)
    print("将改成")
    print("=" * 92)
    _show("after", na)
    print("      提交载荷的键：%s" % sorted(after.keys()))
    if dry:
        print()
        print("DRY-RUN：未提交。加 --go 才真正提交。")
        return
    st, draft = z.req("POST", "/api/records/%s/draft" % rid)
    if st not in (200, 201):
        sys.exit("ERROR: 建 edit 草稿失败 (%s): %s" % (st, json.dumps(draft)[:400]))
    print()
    print("  edit 草稿已建 id=%s" % draft.get("id"))
    payload = {"access": {"record": "public", "files": "public"},
               "files": {"enabled": True},
               "metadata": after}
    st, r = z.req("PUT", "/api/records/%s/draft" % rid, payload)
    if st not in (200, 201):
        sys.exit("ERROR: 写元数据失败 (%s): %s" % (st, json.dumps(r)[:600]))
    print("  元数据已写入 edit 草稿")
    st, r = z.req("POST", "/api/records/%s/draft/actions/publish" % rid)
    if st not in (200, 201, 202):
        sys.exit("ERROR: 发布 edit 失败 (%s): %s" % (st, json.dumps(r)[:600]))
    st, chk = z.req("GET", "/api/records/%s" % rid)
    nc = _norm_meta((chk or {}).get("metadata") or {})
    nc["_rid"] = rid
    nc["title"] = "%s   (doi=%s)" % (nc["title"], (chk or {}).get("doi"))
    print()
    print("  ✅ 订正完成")
    _show("after", nc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--token", default=os.environ.get("ZENODO_TOKEN", ""))
    ap.add_argument("--sandbox", action="store_true")
    ap.add_argument("--net", choices=["auto", "direct", "pin"], default="auto",
                    help="auto=先普通连接，DNS 失败自动切 IP 直连")
    ap.add_argument("--stage", action="store_true", help="建草稿+上传+写元数据（不发布）")
    ap.add_argument("--publish", action="store_true")
    ap.add_argument("--discard", action="store_true")
    ap.add_argument("--edit-record", default="", metavar="ID")
    ap.add_argument("--go", action="store_true", help="--edit-record 时真正提交")
    ap.add_argument("--only", nargs="*", default=[], metavar="KEY",
                    help="只处理这些沉积：raw / apr")
    ap.add_argument("--use", nargs="*", default=[], metavar="ID")
    ap.add_argument("--verify", nargs="*", default=[], metavar="ID",
                    help="核对已 stage 的草稿")
    ap.add_argument("--proxy", default=os.environ.get("ZENODO_PROXY", ""),
                    metavar="HOST:PORT",
                    help="经本地 HTTP 代理的 CONNECT 隧道（实测远快于 pin，"
                         "如 127.0.0.1:5780）")
    args = ap.parse_args()

    if not args.token:
        sys.exit("ERROR: 没给 token。到 https://zenodo.org/account/settings/applications/"
                 "tokens/new/ 建一个（scope 勾 deposit:write），然后\n"
                 "       export ZENODO_TOKEN=...   或   --token ...")
    host = SANDBOX_HOST if args.sandbox else PROD_HOST
    print("=== Zenodo (%s) ===" % ("SANDBOX" if args.sandbox else "PRODUCTION"))
    print("  token : %s" % mask(args.token))
    proxy = None
    if args.proxy:
        ph, _, pp = args.proxy.partition(":")
        proxy = (ph, int(pp or 80))
    print("  host  : %s   net=%s   proxy=%s"
          % (host, args.net, ("%s:%d" % proxy) if proxy else "none"))
    print()
    z = Client(host, args.token, mode=args.net, proxy=proxy)
    st, _ = z.req("GET", "/api/user/records?size=1")
    if st != 200:
        sys.exit("ERROR: token 校验失败 (%s)。检查 scope 是否含 deposit:write。" % st)
    print("  token 有效（实际链路：%s）\n" % (z.chosen or "direct"))

    today = datetime.date.today().isoformat()
    if args.edit_record:
        fix_record(z, args.edit_record, dry=not args.go)
        return
    if args.verify:
        # 只传 id；沉积身份由 verify() 按草稿里的文件名判定，不做位置配对。
        verify(z, [(None, rid, None, None) for rid in args.verify])
        return
    if args.discard:
        if not args.use:
            sys.exit("ERROR: --discard 需要 --use <id> ...")
        discard(z, args.use)
        return
    if args.publish:
        if not args.use:
            sys.exit("ERROR: --publish 需要 --use <id> ...")
        publish(z, args.use)
        return
    if args.stage:
        staged = stage(z, today, only=args.only or None)
        verify(z, staged)
        print()
        print("下一步：浏览器打开上面的草稿页面确认无误，然后")
        print("  python zenodo_deposit.py --publish --use %s"
              % " ".join(str(r) for _, r, _, _ in staged))
        return
    ap.print_help()


if __name__ == "__main__":
    main()
