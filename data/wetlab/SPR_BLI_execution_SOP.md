# 全长 VEGF165 双受体 SPR / BLI 实验执行方案（SOP）
## 抗 VEGF-A 微蛋白候选 1✱/2✱/4✱ 结合与双受体阻断验证（含 #3/#5 探针 + NEG 阴性对照）

> 本 SOP 是 `experimental_design.md` 中 §4（结合实验）的执行版，给出**可直接上机的步骤、缓冲液配方、浓度梯度、动力学拟合与判定矩阵**。
> 所有结合/功能实验**必须使用全长 VEGF165（1–165，含 C 端肝素结合区 HBD）**，严禁用模型截短体（见 §1 红线）。
> 候选蛋白来源：见同目录 `cro_order/CRO_ORDER.md`（His6-SUMO-binder，SUMO 切除后得天然 N 端）。

---

## 1. 中心目标与红线

**三个要回答的问题**
1. 候选是否结合 VEGF165？（是/否）
2. 结合多强？（KD、kon、koff）
3. 是否阻断 VEGF–受体相互作用——且对 **KDR（VEGFR2）** 与 **Flt-1（VEGFR1）** 分别如何？（双受体阻断矩阵）

**红线（不可违反）**
- 🔴 固定相/分析物中的 VEGF 必须是 **全长 VEGF165（His-tagged 重组人，含 HBD 111–165）**。模型仅含 1–98，缺失 HBD 会漏掉 HBD 介导的 KDR D1 / αvβ3 信号（Takada 2024）。
- 🔴 双受体阻断实验必须**独立固定 KDR D2D3 与 Flt-1 D2**，不能用单一受体推断。
- 🔴 计算 ΔG 仅为相对富集值；**最终以 SPR/BLI 的 KD 为准**（experimental_design.md §7.2）。

---

## 2. 物料与缓冲液配方

### 2.1 关键试剂
| 物料 | 规格/来源 | 用途 |
|---|---|---|
| VEGF165（全长，His-tag，人） | Sino Biological 100–32H8 / 自制 | BLI/SPR 固定相、竞争分析物 |
| VEGF110（1–110，无 HBD） | 对照 | 验证 HBD 依赖性 |
| KDR D2D3（Fc 或 His 标签，人） | Sino Biological 100–10H | 双受体阻断（受体芯片） |
| Flt-1 D2（Fc 或 His 标签，人） | Sino Biological 101–08H | 双受体阻断（受体芯片） |
| bevacizumab / aflibercept | 市售 | 阳性对照 |
| scramble 微蛋白 | 同长度同组成重排（本室提供） | 阴性对照 |
| 候选 binder（#1–#5） | 见 CRO_ORDER.md | 分析物 |

### 2.2 缓冲液（实验当天新鲜配制，0.22 µm 过滤）
| 缓冲液 | 配方 | 用途 |
|---|---|---|
| HBS-EP+ | 10 mM HEPES pH 7.4, 150 mM NaCl, 3 mM EDTA, 0.05% v/v P20 | **SPR 运行缓冲液**（Biacore 标准） |
| BLI 运行缓冲 | 1× PBS pH 7.4, 0.1% BSA, 0.02% NaN₃（或 0.01% Tween-20 替代 BSA 视蛋白而定） | **Octet 运行缓冲** |
| 胺偶联试剂 | EDC/NHS（Biacore Amine Coupling Kit）、1 M 乙醇胺 pH 8.5 | CM5 芯片固定 VEGF165 |
| 再生液 A | 10 mM Glycine pH 2.0 | 受体/抗原再生（短接触） |
| 再生液 B | 2 M MgCl₂ | 备选再生（温和） |
| 封闭液 | 1 M 乙醇胺 / 1 mg/mL BSA | 阻断剩余活化位点 |

---

## 3. 阶段 0 — 候选蛋白验收（来自 CRO_ORDER）

上机前必须满足（任一不过则先解决表达/纯化，不进结合实验）：
- [ ] SEC-HPLC 单体峰（无聚集）
- [ ] LC-MS MW 与理论值一致（±0.1%）
- [ ] SUMO 切除确认（binder 天然 N 端）
- [ ] 浓度经 A280 或 AAA 准确定量（误差 < 15%）
- [ ] 非还原 SDS-PAGE + MS 判定二硫键状态（1✱ #1_T0.3_s132 含 2 Cys 为重点；4✱ T0.1_s35 为 0 Cys）

---

## 4. 阶段 1 — CD 折叠验证（1 天，可选但推荐先跑）

- 缓冲：10 mM NaPi pH 7.4 + 20 mM NaCl（低紫外）。
- 样品 0.1–0.5 mg/mL，0.1 mm 光程，25 °C，190–260 nm，1 nm 步长。
- 热变性：222 nm，20→90 °C（2 °C/min）得 Tm。
- **判读**：远紫外有特征谱（α-螺旋 ~208/222 双负峰；β ~218 单负峰）→ 已折叠；随机卷曲 → 重新折叠/复性。Tm > 40 °C 为可开发性底线。

---

## 5. 阶段 2 — BLI 初筛（Octet，高通量，2–3 天）

**策略**：固定 VEGF165（全长）于 AHC/SSA 生物层，浸入 binder 梯度，实时拟合 1:1 动力学。

### 5.1 传感器与固定
- 传感器：**AHC（ streptavidin 型，若 VEGF165 生物素化）** 或 **SSA（如需 His-tag 捕获）**。
- 固定：将 VEGF165（全长，生物素化或 His 标签）上样至 0.5–1 µg/mm² 负载（目标 ~0.5–1 nM 级表面密度，避免 mass transport 限制）。
- 基线：BLI 运行缓冲中平衡 60 s。

### 5.2 结合-解离梯度（每候选）
- **Pilot 梯度（先摸底）**：binder 浓度 0 / 1 / 3 / 10 / 30 / 100 / 300 / 1000 nM（8 点，~3 倍稀释），覆盖假设施 KD 范围。
- 若 Pilot 未饱和 → 加高浓度或更换稀释中心；若过饱和 → 下调起点。
- **方法时序**：基线 60 s → 结合 120–300 s → 解离 300–600 s（解离越长 koff 越准）。
- 每板含：空白传感器（仅缓冲）背景通道 + scramble 阴性 + bevacizumab 阳性（作为动力学质量标杆）。

### 5.3 数据分析（Octet 软件 / Python）
- 扣除背景（参考通道、溶剂校正）。
- 拟合 **1:1 结合模型**（或 1:1 异源二聚/二价若明显偏移）：`R(t) = (kon·C·Rmax)/(kon·C + koff) · (1 − e^{−(kon·C+koff)t})` 形式。
- 输出：**KD = koff/kon**、Rmax、χ²。
- 判定门槛：KD < 1 µM（理想 < 100 nM）进入 SPR 复核；KD > 10 µM 标记"弱/需重设计"。

---

## 6. 阶段 3 — SPR 复核与双受体阻断（Biacore，核心，5–7 天）

### 6.1 通道布置（CM5 芯片，胺偶联固定 VEGF165）
| 通道 | 固定相 | 用途 |
|---|---|---|
| FC1 | VEGF165（全长） | 直接结合（分析物 = binder 梯度） |
| FC2 | KDR D2D3 | 双受体阻断（竞争） |
| FC3 | Flt-1 D2 | 双受体阻断（竞争） |
| FC4 | BSA / 空白 | 背景/非特异性 |

- 胺偶联：EDC/NHS 活化 → 上样 VEGF165（~5–10 µg/mL in 10 mM NaAc pH 4.5，~500–1000 RU）→ 乙醇胺封闭。
- 受体通道（FC2/FC3）：固定 KDR D2D3 / Flt-1 D2 至 ~300–500 RU。
- 对照：单独固定 bevacizumab 已知结合 VEGF165，验证芯片活性。

### 6.2 直接结合动力学（FC1，多循环 MCK）
- binder 梯度同上（Pilot 后定标），流速 30–50 µL/min，23 °C。
- 循环：结合 60–120 s → 解离 120–300 s → 再生（Gly pH 2.0，≤30 s 接触，勿伤 VEGF165）。
- 拟合 1:1 / 双态模型，得 KD/kon/koff；与 BLI 交叉验证（偏差 > 3× 需排查 mass transport 或固定密度）。

### 6.3 双受体阻断实验（关键，对应结合模式结论）
**竞争法**（每候选）：
1. 预孵育：VEGF165（固定浓度，如 10 nM）+ binder（梯度 0 / 0.1 / 1 / 10 / 100 × IC₅₀ 估）于运行缓冲，室温 30 min。
2. 将预孵育混合物流经 **FC2（KDR）** 与 **FC3（Flt-1）**，记录受体结合信号。
3. 计算**阻断率** = (1 − R_sample/R_VEGF_only) × 100%。
4. 以 binder 浓度对阻断率拟合，得 **IC₅₀(阻断)**。

**判定矩阵（结合模式预测 → 实验验证）**：

| 候选 | 预测阻断 KDR? | 预测阻断 Flt-1? | 实验判定标准 |
|---|---|---|---|
| 4✱ (T0.1_s35) | **是**（命中 KDR 位点） | 可能部分 | KDR 通道阻断 IC₅₀ < binder KD 量级 → 确认主靶点；Flt-1 通道若弱阻断则选择性成立 |
| 1✱ (#1_T0.3_s132), #3 (Flt-1 探针) | 待验证 | **是** | **必须证明确实交叉阻断 KDR**；若 KDR 通道零阻断 → 效价存疑，降级 |
| 2✱ (#2c_T0.1_s18, 双受体) | **是**（双受体 Flt-1/KDR 区） | 可能部分 | 期望 KDR + Flt-1 双通道均显著阻断（广谱）→ 本候选最大卖点 |
| #5 (vegf_len100_3, Flt-1 主导) | 是（Flt-1 主导） | 是 | Flt-1 通道强阻断、KDR 通道**零/弱阻断**（KDR 0 残基）→ 定位为 Flt-1 选择性候选，不阻断主效价信号 |

> 阳性对照 bevacizumab/aflibercept 应双通道显著阻断；阴性 scramble 应无阻断。

---

## 7. 阶段 4 — 细胞功能验证（与 SPR 互证，1–2 周）

- **HUVEC 增殖**：VEGF165（10–20 ng/mL）+ binder 梯度，CCK-8 / EdU，求 IC₅₀。
- **管形成**：Matrigel + HUVEC + VEGF165 ± binder。
- **受体磷酸化 Western**：p-VEGFR2(Tyr1175) / p-VEGFR1(Tyr1333) —— **直接判定阻断 KDR 还是 Flt-1 信号**，与 §6.3 互为印证（SPR 阻断 KDR ⇔ p-VEGFR2↓）。
- 报告细胞：稳定表达 KDR 或 Flt-1 工程细胞（如有）。

---

## 8. 数据分析与拟合要点

- 至少一个独立重复（n≥2），报告均值 ± SD。
- 动力学用 **1:1 Langmuir** 为首选；若 sensorgram 呈明显 hook / 双相，试 **两态（conformational）** 或 **异源二聚（binder 二价于 VEGF 二聚体）** 模型并比较 AIC。
- 竞争阻断用 **4-参数 logit** 拟合阻断率–浓度曲线。
- 全部 KD 与 SPR/BLI 阻断 IC₅₀ 汇总入 `results/candidates_6_final.md` 实验栏。

---

## 9. Go / No-Go 决策门槛

推进 lead optimization 须同时满足：
- SPR/BLI **KD < 1 µM**（理想 < 100 nM）；
- **至少阻断 KDR**（或 KDR+Flt-1 双阻断）；#5 为 Flt-1 选择性特例（须明确标注不阻断 KDR 信号）；
- 细胞 **IC₅₀ < 10 µM** 且 p-VEGFR 抑制明确；
- 单体、可溶、Tm > 40 °C。

**优先级执行顺序**（资源约束时）：
**4✱ 先行**（KDR 位点 T0.1_s35，药效学最相关）→ **1✱ / 2✱ 并行**（综合分 Top-2，须过双受体关）→ **#3 跟进** → **#5 备选**（Flt-1 选择性专项）。全程含 NEG 阴性对照。

---

## 10. 时间线与排程

| 周次 | 阶段 | 交付 |
|---|---|---|
| W1 | 阶段 0 验收 + 阶段 1 CD | 折叠确认、Tm、可上机判定 |
| W2–W3 | 阶段 2 BLI 初筛（5 候选 + 对照） | 各候选 KD（初值）、结合/非结合 |
| W3–W4 | 阶段 3 SPR 直接结合 + 双受体阻断 | KD/kon/koff（复核）、KDR/Flt-1 阻断 IC₅₀ |
| W4–W6 | 阶段 4 细胞功能 | IC₅₀、p-VEGFR Western |
| W6 | 决策 | Go/No-Go 报告，确定 lead |

---

## 11. Troubleshooting

| 现象 | 可能原因 | 对策 |
|---|---|---|
| BLI 信号漂移/无结合 | VEGF165 固定失败或失活 | 检查 Biotin/His 标记效率；换批；用 bev 阳性验芯片 |
| SPR 基线噪声高 | 固定密度过高 / 再生不当 | 降固定 RU；缩短再生接触；换 MgCl₂ 再生 |
| 仅解离慢、kon 不稳 | mass transport 限制 | 降表面密度、升流速至 100 µL/min |
| KDR 通道零阻断但预测应阻断 | binder 仅结合 VEGF N 端而受体结合区在别处，或全长构象差异 | 补 VEGF110 vs VEGF165 对照；复核结合模式残基 |
| 候选聚集（SEC 多峰） | 多 Cys 错误配对 / 疏水暴露 | 1✱ (#1_T0.3_s132, 2 Cys) 用氧化折叠菌株；单 Cys 候选考虑 Cys→Ser 去聚集（见 CLONING.md §5） |

---

## 12. 关联文件

- 实验设计总纲：`experimental_design.md`
- 结合模式/表位：`../../../results/binding_mode/binding_mode_report.md`
- CRO 克隆工单：`cro_order/CRO_ORDER.md`
- 构建体/引物：`constructs/CLONING.md`、`constructs/constructs_table.csv`
- 综合终表：`../../../results/candidates_6_final.md`
