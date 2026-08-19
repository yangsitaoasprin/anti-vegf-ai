# Upload guide — GitHub (code) + Zenodo (raw data)

This repository is **local only** so far. Two public deposits are required for CSBJ
compliance (Elsevier does **not** accept "available on request").

> You (the author) must perform the actual push / upload — they need your credentials.
> This guide tells you exactly what to do. After each deposit, paste the resulting DOI
> back into `README.md`, `results_manifest.md`, `CODE_AUTHORS.md`, and the manuscript's
> Data/Code availability statement.

---

## 1. GitHub — the code & minimal data in this `repro/` folder (step by step)

The `repro/` folder is small (~1.4 MB of bundled data + scripts) — well within GitHub
limits. **Do not** `git add` the 311 MB `results/` raw outputs here (those go to Zenodo).

> The local `repro/` is already `git init`-ed and committed (52 files). You only need to
> **create the remote repo and push**.

### 1.1 准备（已做，确认即可）
```bash
cd /path/to/anti_VEGF_AI_project/repro
git log --oneline -1        # 应看到 "Reproducible code & minimal data ..." 提交
git status                  # 应显示 "nothing to commit, working tree clean"
```

### 1.2 在 GitHub 上建空仓库（三种方式任选）

**A. 网页（最直观）**
1. 打开 https://github.com/new
2. Repository name 建议：`anti-vegf-ai`
3. **Visibility 选 Public**
4. **不要**勾 "Add a README" / ".gitignore" / "license"（本地已有，避免冲突）
5. 点 **Create repository**

**B. gh CLI（若已装并登录）**
```bash
gh repo create anti-vegf-ai --public --description "Reproducible code & minimal data for anti-VEGF-A miniprotein virtual screen" --source . --remote origin --push
# 上面一行会直接建仓并 push；若只想建仓不 push，去掉 --push
```

**C. GitHub Desktop（图形界面，免命令行）**
1. 打开 GitHub Desktop → **File → Add Local Repository** → 选 `repro/` 目录
2. **Publish repository** → Name 填 `anti-vegf-ai`、勾 Public → **Publish**
3. 之后本地改动用 **Push origin** 按钮同步

### 1.3 关联远程并推送（方式 A/B 的命令行补完）
```bash
git branch -M main
git remote add origin https://github.com/yangsitaoasprin/anti-vegf-ai.git
git push -u origin main
```
> 若提示认证：GitHub 自 2021 起**不能用账号密码**，用 **Personal Access Token (PAT)**
> 代替密码（Token 需在 https://github.com/settings/tokens 生成，勾 `repo` 权限）；
> 或配置 SSH key（`git@github.com:...` 形式 remote）更省事。

### 1.4 拿 GitHub DOI（软件引用用）

CSBJ 接受两种形式，任选其一写进声明：

- **方式一（推荐，得到正式 DOI）**：通过 **Zenodo 的 GitHub 集成** 给代码仓库发 DOI。
  > ⚠️ GitHub 自身**不发 DOI**。下面拿到的 DOI 实际由 **Zenodo** 签发，只是它绑定到你的 GitHub release 快照。
  > 这也是回填脚本里 `<github-doi>` 的真实来源——它和原始数据 DOI 一样都是 `10.5281/zenodo.XXXXXXX` 格式。

  1. 先在 Zenodo 侧授权仓库（只需做一次）：
     登录 https://zenodo.org → 右上角头像 → **GitHub**（或直接开
     `https://zenodo.org/account/settings/github/`）→ 点 **Authorize Zenodo / Connect**
     用 GitHub 账号授权。
     - 授权后页面列出你 GitHub 账号下**所有仓库**，每个仓库一行，行尾有一个 **开关 toggle**
       （默认关）。把 `anti-vegf-ai` 的开关**打开**——开启后 Zenodo 会监听该仓库的 **release 事件**。
     - 开关右侧还有一个 **Sync / Update** 按钮：当仓库被改名、transfer 或刚 push 了新 release 但
       Zenodo 没抓到时，点它手动重新同步一次。
     - 授权成功后 GitHub 仓库主页会出现一个 **Zenodo** 徽章（或在 Settings → Integrations 里看到
       Zenodo 已连接）。
     - ⚠️ **仓库必须是 public**：Zenodo 无法抓取 private 仓库；若仓库设了私有，请先在 GitHub
       设为 Public（Settings → Change visibility）再开启开关。
     - ℹ️ 此时**还没有 DOI**——DOI 要等你打了 release（下一步）后 Zenodo 才生成。
  2. 在 GitHub 给仓库打一个 **tagged release**：
     ```bash
     git tag v1.0
     git push origin v1.0
     ```
     在 GitHub 页面 **Releases → Draft a new release**，选 `v1.0`，填标题/说明，点 **Publish release**。
  3. Zenodo 会自动检测到这个新 release 并抓取快照，生成一条新 deposit；
     稍等片刻（通常几分钟内）刷新 **Zenodo → Upload → 我的上传**，找到对应记录。
  4. 进该记录页，点 **Publish**（首次可能需补一行标题/作者），发布后顶部出现
     `10.5281/zenodo.XXXXXXX` —— **这就是你的 GitHub DOI**，复制它（+ `https://doi.org/...`）。

- **方式二**：直接用仓库 URL（`https://github.com/yangsitaoasprin/anti-vegf-ai`）作为代码位置，
  在声明里注明 "code repository (accessed <日期>)"。CSBJ 接受这种，但带 DOI 更规范。

### 1.5 回填
把仓库 URL（和/或 GitHub DOI）填回第 3 节列出的占位处。

---

## 2. Zenodo — the 311 MB raw model outputs (step by step)

Zenodo (https://zenodo.org) is a free, CERN-operated data repository that mints a
permanent **DOI** for any upload. It is the right home for the 311 MB raw model
outputs that GitHub cannot hold. Below is the exact click-by-click flow.

### 2.1 准备上传包（二选一）

**A. 直接上传目录（推荐，Zenodo 会自动打成 zip）**
只需把 `results_manifest.md` 里列出的目录准备好。不要上传日志/废弃文件。

**B. 先压缩成一个 tar.gz（更可控、断点友好）**
```bash
cd /path/to/anti_VEGF_AI_project
# 仅打包论文真正依赖的原始输出（见 results_manifest.md 的目录清单）
tar czf anti_vegf_raw_results.tar.gz \
    results/deA_redesign_1_2 \
    results/mmgbsa_redesign \
    results/candidates_6_final.csv
# 校验体积
ls -lh anti_vegf_raw_results.tar.gz
```
> 单文件上限 50 GB，311 MB 完全没问题。若担心网络，可拆成 2–3 个 tar 分卷。

### 2.2 注册 / 登录

1. 打开 https://zenodo.org → 右上角 **Log in**。
2. 推荐选 **"Log in with GitHub"**（用你上传代码的同一个 GitHub 账号，便于关联）。
3. 首次登录会要求填 **ORCID**（可选，建议填，对学术履历有利）和确认邮箱。

### 2.3 新建 deposit 并上传文件

1. 顶部 **Upload** → **New upload**（或右上角 **+** → **New upload**）。
2. 进入编辑页后，直接把 `anti_vegf_raw_results.tar.gz`（或整个 `results/` 目录）**拖入**虚线框，
   或点 **Choose files**。大文件会自动分块上传，进度条走完前**不要关闭页面**。
3. 可顺手再拖一个 `results_manifest.md` 进去，作为文件清单说明。

### 2.4 填写元数据（每一项都对照填）

| 字段 | 填什么 |
|---|---|
| **Upload type** | `Dataset` |
| **Publication date** | 保持默认（今天）或选投稿日 |
| **Title** | *Raw model outputs for "Developability-gated, multi-model virtual screening of de novo anti-VEGF-A miniprotein binders"* |
| **Authors** | `Yang Sitao`； affiliation 填 `School of Pharmaceutical Sciences, Dali University; The Third Affiliated Hospital of Dali University`；可加 ORCID |
| **Description** | 复制 `results_manifest.md` 的摘要段：说明包含 Boltz-1 pdb、OpenFold3 cif/json、MM-GBSA 逐种子值，及与论文 Table 2 / Fig. 4 的对应关系 |
| **License** | **Creative Commons Attribution 4.0 (CC-BY-4.0)**（CSBJ 偏好开放许可；不要用 "Restricted") |
| **Keywords** | `VEGF-A`, `miniprotein binder`, `Boltz-1`, `OpenFold3`, `MM-GBSA`, `de novo protein design` |
| **Related identifiers** | 点 **+ Add**，填 GitHub 仓库 URL，`Relation` 选 **"Is supplemented by"**（或 "References"），`Resource type` 选 `Software` |
| **Communities** | 可选：加入 `biosciences` 或你单位的 Zenodo community（非必须） |
| **Grants / Funding** | 若论文 M4 填了基金号，在此关联（非必须） |

> 勾选 **"Reserve DOI"** 之前的预览页也能看到，但 DOI 只在 **Publish** 后才正式生效。

### 2.5 发布并拿到 DOI

1. 确认文件 + 元数据无误后，点右下角 **Publish**。
2. 弹窗确认 → 再次 **Publish**。
3. 发布成功后页面顶部出现 **`10.5281/zenodo.XXXXXXX`**（XXXXXXX 是数字，例如 `10.5281/zenodo.1423456`）。
4. **复制这个 DOI**（连同完整链接 `https://doi.org/10.5281/zenodo.XXXXXXX`）。

### 2.6 之后要更新数据怎么办（版本管理）

- Zenodo 的 DOI 是**版本化**的：点 **New version** 可上传修订，旧 DOI 仍指向旧版，新版本获得新 DOI，
  同时有一个**总 DOI**（`10.5281/zenodo.XXXXXXX` 无后缀）始终指向最新版。
- 论文里建议引用**总 DOI（无版本号后缀）**，这样将来你补数据，读者永远拿到最新版。

### 2.7 回填占位

把 `10.5281/zenodo.XXXXXXX` 填回下面几处（见第 3 节）：
`README.md`、`results_manifest.md`、`CODE_AUTHORS.md`、以及两版稿件的 Data/Code availability。

---

## 3. Backfill the DOIs

Replace the placeholders (do this **after** you have both DOIs from §1.4 and §2.5):

| File | Placeholder to replace |
|---|---|
| `README.md` | Zenodo DOI line |
| `results_manifest.md` | `10.5281/zenodo.XXXXXXX` |
| `CODE_AUTHORS.md` | `<zenodo-DOI>`, `<repo-URL>` |
| `论文/Manuscript_v7.md` (Data/Code availability) | GitHub URL + Zenodo DOI |
| `论文/Manuscript_ZH.md` (数据/代码可用性) | 同上（中文） |

> **IMPORTANT — the two `XXXXXXX` mean different things.** In the two manuscripts,
> `10.5281/zenodo.XXXXXXX` is the **GitHub code-repo DOI** (from §1.4); in
> `results_manifest.md` it is the **Zenodo raw-data DOI** (from §2.5). Do not paste the
> same value into both.

### 3.1 One-shot helper (recommended)

`backfill_dois.py` in this folder does all replacements at once, scoped per file so the
two `XXXXXXX` never get mixed up. It is **dry-run by default** (prints what would change,
writes nothing) — add `--apply` to actually write:

```bash
cd /path/to/anti_VEGF_AI_project/repro

# preview (safe, no writes)
python backfill_dois.py \
    --repo-url https://github.com/yangsitaoasprin/anti-vegf-ai \
    --github-doi 10.5281/zenodo.<代码DOI号> \
    --zenodo-doi 10.5281/zenodo.<原始数据DOI号>

# apply the changes
python backfill_dois.py ... --apply

# commit the backfilled files, then push
git add README.md results_manifest.md CODE_AUTHORS.md \
         ../论文/Manuscript_v7.md ../论文/Manuscript_ZH.md
git commit -m "Backfill GitHub + Zenodo DOIs"
git push
```

The script prints a `git diff --stat` after `--apply` so you can verify the edits.
(Paths are relative to `anti_VEGF_AI_project/`, so the repo works wherever it is cloned.)

### 3.2 Manual edit (alternative)

If you prefer to edit by hand, open each file in the table above and replace the
placeholder strings with the corresponding values from §1.4 / §2.5.

---

## 4. CSBJ Data/Code availability — required wording (already drafted in the manuscript)

> *Data availability.* Representative complex structures, the template (PDB 1FLT), the
> final candidate table and the per-seed MM-GBSA decomposition are available in the
> GitHub repository at `<repo-URL>` (DOI: `<github-DOI>`). The complete raw Boltz-1 and
> OpenFold3 predictions (~311 MB) are deposited on Zenodo at
> `https://doi.org/<zenodo-DOI>` under CC-BY-4.0. VEGF-A structural data derive from PDB
> entry 1FLT. Code to regenerate all figures and tables is in the same GitHub repository
> (MIT license).
