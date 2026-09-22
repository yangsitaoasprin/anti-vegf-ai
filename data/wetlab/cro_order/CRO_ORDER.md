# CRO 克隆工单 — 抗 VEGF-A 微蛋白候选 #1–#5（His6-SUMO / pET-28a(+)）

> 配套序列/引物文件见同目录：
> - `fusion_protein.fasta` — His6-SUMO-binder 融合蛋白氨基酸序列
> - `insert_orf_dna.fasta` — 密码子优化插入 ORF（NdeI…XhoI，含起始/终止）
> - `constructs_cro.csv` — 结构化清单（含融合序列/ORF DNA/引物/二硫键位点）
> - `primers_cro.fasta` — 限制酶（NdeI/XhoI）+ Gibson 双版引物
>
> 生成日期：2026-07-28｜关联设计：`../experimental_design.md`、结合模式 `../../../results/binding_mode/binding_mode_report.md`

---

## 1. 工单概要

| 项目 | 内容 |
|---|---|
| 目标 | 获得 5 个候选（#1–#5）可表达、可纯化的重组融合蛋白：N 端 His6–SUMO–binder |
| 表达宿主 | *E. coli* BL21(DE3)（CRO 或本室转化） |
| 载体 | **pET-28a(+)**（Novagen/Merck），NdeI / XhoI 多克隆位点，C 端 His 标签回路 |
| 标签策略 | N 端 **His6–Smt3(酵母 SUMO)**，SUMO 蛋白酶（Ulpl1/SENP2）切去，得**天然 N 端 binder** |
| 交付形式 | ① 全基因合成 His6-SUMO-binder 插入 ORF（密码子优化，E. coli）；② 克隆入 pET-28a(+) 的重组质粒（测序验证）；或 ③ 仅合成基因+引物由本室自克隆 |
| 优先级 | 4✱（KDR 位点 T0.1_s35，最高）→ 1✱/2✱（综合分 Top-2，双受体验证）→ #3 → #5（Flt-1 选择性） |

---

## 2. 载体与标签策略（关键纠错说明）

**本工单使用的 SUMO 为正确的酵母 Smt3（UniProt Q12306）**，成熟序列 1–98 以 **Gly-Gly（GG）** 结尾，可被 SUMO 蛋白酶精确切除，切后剩余 binder 天然 N 端（无多余残基）：

```
Smt3 (酵母, 1-98): SDSEVNQEAKPEVKPEVKPETHINLKVSDGSSEIFFKIKKTTPLRRLMEAFAKRQGKEMDSLRFLYDGIRIQADQTPEDLDMEDNDIIEAHREQIGG
                                                                                                         └─ GG (蛋白酶切点)
```

- **His6 标签**：`MHHHHHH`（起始 M 由 NdeI 的 ATG 提供，故融合构建体写作 `MHHHHHH` + `SDSEVN...`）。
- **融合构型**：`Met – His6(6×) – Smt3(98aa) – binder(60–100aa)`，全长 164–204 aa。
- ⚠️ **避免使用原 `design_constructs.py` 中的 SMT3 变量**：该序列 C 端非 Gly-Gly（蛋白酶无法切割），**以本工单为准**。
- ⚠️ **不使用 pE-SUMO（LifeSensors）的 BsaI Golden-Gate 方案**：其连接区为专有序列，公开渠道无法获得确切同源臂；改用 pET-28a(+) 的**公开确定同源臂**（NdeI/XhoI），更稳妥可投递。

---

## 3. pET-28a(+) 真实同源臂（公开序列，用于克隆/引物设计）

| 位置 | 同源臂序列（5′→3′） | 说明 |
|---|---|---|
| 5′ 上游（T7 RBS → NdeI 起始） | `GAAATAATTTTGTTTAACTTTAAGAAGGAGATATACATATG` | 末位 `CATATG` = **NdeI**（ATG 为起始 M） |
| 3′ 下游（XhoI → 载体 C-His/终止区） | `CTCGAGCACCACCACCACCACCACTGAGATCCGGCTGCTAACAAAGCCCGAAAGG` | 前 6 位 `CTCGAG` = **XhoI** |

- 插入 ORF 构型（合成基因）：`CATATG` + `ATG`+`His6`+`SUMO`+`binder`（密码子优化）+ `TAA` + `CTCGAG`
- 克隆后表达框：`T7 promoter – RBS – CATATG(ATG-His6-SUMO-binder)–TAA–CTCGAG – 终止子`
- binder C 端后接载体自带 **C 端 His6**（由 XhoI 下游回路表达），便于 Ni-NTA 捕获；SUMO 切除后如需去除 C-His 需另设 TEV/蛋白酶位点（本工单默认保留 C-His 用于纯化，不影响功能验证）。

---

## 4. 候选清单（综合分终表排名）

| 排名 | 候选 ID | binder 长度 | 融合总长 | 命中 VEGF 表位 | Cys 数(位置) | 克隆方式 |
|---|---|---|---|---|---|---|
| #1 | #1_T0.3_s132 (1✱) | 80 | 184 aa | Flt-1 + 二聚面 | 2 (7,75) | NdeI/XhoI 或 Gibson |
| #2 | #2c_T0.1_s18 (2✱) | 100 | 204 aa | KDR + Flt-1 + 二聚面（双受体广谱） | 2 (18,86) | NdeI/XhoI 或 Gibson |
| #3 | vegf_len60_1_T0.2_s1 | 60 | 164 aa | Flt-1 + 二聚面（机制探针） | 2 (19,38) | NdeI/XhoI 或 Gibson |
| #4 | T0.1_s35 (4✱) | 80 | 184 aa | **KDR（药效学最相关，最高优先）** | 0 | NdeI/XhoI 或 Gibson |
| #5 | vegf_len100_3_T0.1_s2 | 100 | 204 aa | Flt-1 主导（KDR 0，选择性探针） | 1 (92) | NdeI/XhoI 或 Gibson |
| #0 | NEG_scrambled_T0.1_s35 | 80 | 184 aa | 阴性对照（序列打乱） | 0 | NdeI/XhoI 或 Gibson |

> **克隆方式说明**：全部 6 个构建体（含 #0 NEG）ORF 内部**均无 NdeI/XhoI 位点**，故 NdeI/XhoI 限制酶克隆与 Gibson 组装两种方式均可，CRO 可任选。早期版本中「#4 含内部 NdeI 须走 Gibson」是针对旧构造 vegf_len80_3 的注记，现已由去淀粉样替代版 T0.1_s35（0 Cys、无内部酶切位点）取代，该限制不再适用。

---

## 5. 引物表（与 constructs_cro.csv 完全一致，NdeI/XhoI + Gibson）

### 5.1 限制酶克隆引物（NdeI / XhoI）

| 候选 | 正向（NdeI_FWD，共用） | 反向（XhoI_REV） |
|---|---|---|
| 全部 | `AAGAAGGAGATATACATATGCATCATCATCATCATCATAG`（RBS + NdeI CATATG + His6 起始） | — |
| #1 | — | `CCGCTCGAGTTACGCTTCCATCGCCTGGCATTC` |
| #2 | — | `CCGCTCGAGTTACAGCGCCGCCGCCTGCGCTTT` |
| #3 | — | `CCGCTCGAGTTAACGCAGCTGACGCACCAGTTC` |
| #4 | — | `CCGCTCGAGTTAACGTTTCGCCGCCACCGCCGC` |
| #5 | — | `CCGCTCGAGTTATTTCGCCAGTTTCGCCGCTTC` |

### 5.2 Gibson 组装引物（pET-28a NdeI/XhoI 线性化）

| 候选 | 正向（GIBSON_FWD，共用） | 反向（GIBSON_REV） |
|---|---|---|
| 全部 | `GGGAATTCCATATGCATCATCATCATCATCATAGC`（EcoRI + NdeI + His6 起始） | — |
| #1 | — | `CGCTTCCATCGCCTGGCATTTGGTGGTGGTGGTGCTCGAG` |
| #2 | — | `CAGCGCCGCCGCCTGCGCTTTGGTGGTGGTGGTGCTCGAG` |
| #3 | — | `ACGCAGCTGACGCACCAGTTTGGTGGTGGTGGTGCTCGAG` |
| #4 | — | `ACGTTTCGCCGCCACCGCCGTGGTGGTGGTGGTGCTCGAG` |
| #5 | — | `TTTCGCCAGTTTCGCCGCTTTGGTGGTGGTGGTGCTCGAG` |

> 完整序列（含融合蛋白 AA、ORF DNA）见 `fusion_protein.fasta` / `insert_orf_dna.fasta` / `constructs_cro.csv`。


## 6. 克隆执行协议（供 CRO 或本室）

1. **基因合成**：合成 6 条插入 ORF（#1–#5 + #0 NEG，见 `insert_orf_dna.fasta`），5′ 端 `CATATG`（含 ATG）、3′ 端 `TAA`+`CTCGAG`，E. coli 密码子优化（全部已避开内部 NdeI/XhoI/NcoI/BamHI 禁忌）。
2. **载体线性化**：pET-28a(+) 用 **NdeI + XhoI** 双酶切（全部 6 个候选 ORF 内部均无 NdeI/XhoI，两种方式均可）。
3. **连接/组装**：
   - 限制酶路径：T4 DNA 连接酶，16 °C 过夜；全部候选。
   - Gibson 路径：线性化载体 + 插入 ORF，Gibson Assembly Master Mix，50 °C 1 h；全部候选亦可。
4. **转化**：DH5α 克隆株，卡那霉素（50 µg/mL）筛选。
5. **测序验证**：双向 Sanger，覆盖全长插入（重点：NdeI–SUMO–binder 连接区、binder 全长无突变）。
6. **表达株转化**：验证无误质粒转 BL21(DE3)，甘油菌保藏。

---

## 7. 表达 / 纯化 / 酶切 QC（验收标准）

| 步骤 | 条件 | 验收 |
|---|---|---|
| 诱导表达 | 0.5 mM IPTG，18 °C，16–18 h（低温降聚集） | 可溶性组分含 ~17–23 kDa 融合带（His6-SUMO-binder） |
| 捕获 | Ni-NTA，咪唑 20→250 mM 梯度洗脱 | 单一条带 |
| 精纯 | SEC-HPLC（Superdex 75），20 mM NaPi / 150 mM NaCl pH 7.4 | 单体峰（非二聚/聚集） |
| SUMO 切除 | Ulp1 / SENP2（1:50 w/w），4 °C 过夜 | binder 天然 N 端释放，经 Ni-NTA 去除 His6-SUMO 标签 |
| 质谱 | LC-MS / ESI | 实测 MW 与理论值一致（±0.1%） |
| 二硫键 | 非还原 SDS-PAGE + MS | #1（4 Cys）需氧化折叠判定；单 Cys 候选确认非错误配对 |
| 定量 | A280（#2/#4/#5 含 Trp 或按序列）；#1/#3 无 Trp → AAA 或定量肽图 | 浓度误差 < 15% |
| **聚集 QC（#2 与 #4）** | SEC-HPLC 单体率 + 浊度 A340（每批必测） | 单体率达标方可进 SPR；不达标则该批重表达或退回设计。
| 　└ #2 触发理由 | 官方 APR-Score nHigh 12/95，六个构建体最高，高于阴性对照 4/75 | 计数受缺失特征 SH_ENTR 影响（区间 4–20/95），但 #2 在区间内始终最高 |
| 　└ #4 触发理由 | APR 计数不高（3/75），但结构层显示 3 个窗口全部表面暴露（31.0 % / 31.4 % / 35.9 %）、不在结合界面、靶点不遮蔽，且与父本 WALTZ 高风险区同位 | 暴露度为单一构象计算，未做 MD 采样 |

> **#2 专项聚集 QC 的依据（2026-09-14 增补）**：用**官方 APR-Score 权重**（LR+RF 共识，阈值 P>0.5）对全部构建体重算，`#2c_T0.1_s18` 的 nHigh = **12/95**（受缺失特征 SH_ENTR 影响，官方公布区间内为 **4–20/95**；但 #2 在该区间内**始终最高**，构建体之间排序不变），是六个构建体中最高的一个，也是唯一高于阴性对照的候选（NEG 4/75；#1 0/75；#4 3/75；#3 2/55；#5 8/95）。该风险来自 #2 的第二步"界面保守微调"（H72Y + A85E + K95A），不是第一步去丙氨酸（去 A 步骤把 ≥4 Ala 的六肽窗口从 15 降到 0）。据 `results/apr_official_integration/beforeafter_report.txt`。

---

## 8. 对照与配套物料

- **阳性对照**：bevacizumab（Avastin）、aflibercept（Eylea）— VEGF165 结合/功能阳性。
- **阴性对照**：同长度同组成 scramble 微蛋白（本室另提供）验证特异性。
- **全长抗原（关键）**：**VEGF165（重组人，全长 1–165，含 C 端肝素结合区 HBD）** — 所有结合/功能实验必须用全长，不可用模型截短体（见 experimental_design.md §7）。
- **双受体蛋白**：KDR D2D3、Flt-1 D2 胞外域（Sino Biological 或自制）。

---

## 9. 风险与注意

1. **无内部酶切位点禁忌**：全部 6 个候选（含 #0 NEG）ORF 内部均无 NdeI/XhoI，NdeI/XhoI 限制酶克隆与 Gibson 均可，CRO 任选（旧「#4 含内部 NdeI 须走 Gibson」针对已弃用的 vegf_len80_3，不适用）。
2. **二硫键**：模型为还原态几何，真实氧化态需实验判定；1✱ (#1_T0.3_s132) 含 2 Cys 建议氧化折叠菌株（如 Origami/SHuffle）或表达后体外氧化复性；4✱ (T0.1_s35) 为 0 Cys，无需处理。
3. **可开发性未评分**：本工单仅覆盖"获得可纯化蛋白"，T 细胞表位/聚集/蛋白酶稳定性需下游补评。
4. **序列以本工单文件为准**：`constructs_cro.csv` / `*.fasta` 为权威；请勿使用旧 `design_constructs.py` 的 SMT3 序列。

---

## 10. 投递清单（checklist）

- [ ] 6× 插入 ORF DNA（#1–#5 + #0 NEG，`insert_orf_dna.fasta`，已优化）
- [ ] 6× 融合蛋白 AA（`fusion_protein.fasta`）
- [ ] 引物：NdeI_FWD（1 条共用）+ XhoI_REV（5 条）+ GIBSON_FWD（1 条共用）+ GIBSON_REV（5 条）
- [ ] 载体：pET-28a(+)
- [ ] 克隆方式标注：全部 = NdeI/XhoI 或 Gibson（无内部位点禁忌）
- [ ] 验收：测序 + SEC 单体 + LC-MS MW 一致 + SUMO 可切除
