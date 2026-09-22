# 抗 VEGF-A 微蛋白候选 #1–#5 实验设计
## 合成 / 圆二色谱(CD) / SPR·BLI 结合实验

> 配套文件：
> - `candidates_1_5_sequences.fasta` — 合成用序列（含注释头）
> - `candidates_1_5_for_synthesis.csv` — 结构化清单
> - `structures/<id>.complex.pdb` — Boltz-1 预测结合复合物（链A=VEGF，链B=binder）
> - `structures/*_fullVEGF_*_complex.pdb` — **全长 VEGF165 复核复合物**（链A=全长 VEGF165(165aa)，链B=binder），供实验室参考最新预测结构
> - `../results/binding_mode/binding_mode_report.md` — 结合模式与表位分析（决定双受体验证策略）
> - `../results/deA_redesign_1_2/report_fullVEGF_recheck.md` — 全长 VEGF165 Boltz + MM-GBSA 复核报告

---

## 1. 候选面板与优先级（已对齐最终 CRO 投递集 1✱/2✱/4✱）

> **重要**：下表 ID 已与 `cro_order/constructs_cro.csv` 及 `results/candidates_6_final.md` 对齐。
> 早期探索性候选（vegf_len80_1 / vegf_len100_4 / vegf_len80_3 / vegf_len60_1 / vegf_len100_3）经
> **de-A 去聚集重设计（#1→#1_T0.3_s132）**、**保守界面微调（#2→#2c_T0.1_s18）**、**去淀粉样突变（#4→T0.1_s35）**
> 升级，下表为最终送样/实验版。计算验证（含**全长 VEGF165 复核**）已全部通过，见 `results/deA_redesign_1_2/report_fullVEGF_recheck.md`。

| 湿实验编号 | 对应计算候选 | 长度 | 设计定位 | ipTM(98aa) | ΔG(kcal/mol) | 全长复核 ΔG(干净) | 命中 VEGF 表位 | 优先级 |
|---|---|---|---|---|---|---|---|---|
| 1✱ | #1_T0.3_s132 (de-A) | 80 | #1 去A重设计优选 | 0.937 | −129.9 | **−96.8** | Flt-1(VEGFR1)位点 + VEGF二聚面 | **最高(结合最强)** |
| 2✱ | #2c_T0.1_s18 (A85E+K95A) | 100 | 保守界面微调 | 0.628 | −53.3 | **−65.7** | **KDR+Flt-1+二聚面（双受体广谱）** | 高(双受体，综合关键) |
| 4✱ | T0.1_s35 (16突变) | 80 | #4 去淀粉样替代 | 0.699 | −51.4 | **−95.4** | **KDR(VEGFR2)位点** | 高(药效学最相关) |
| #3 | vegf_len60_1_T0.2_s1 | 60 | Flt-1 探针（对照） | 0.739 | −49.9 | — | Flt-1 位点 + 二聚面 | 中(机制探针) |
| #5 | vegf_len100_3_T0.1_s2 | 100 | Flt-1 选择性探针 | 0.603 | −38.3 | — | Flt-1 主导（KDR 0） | 中(选择性探针) |
| NEG | NEG_scrambled_T0.1_s35 | 80 | 阴性对照（序列打乱） | — | — | — | 无 | 必做(排除假阳) |

**关键结论（结合模式 + 全长复核）**：
- VEGF 主要促血管生成信号由 **KDR/VEGFR2** 介导 → **4✱ (T0.1_s35) 命中 KDR 位点，药效学最相关**。
- **2✱ (#2c_T0.1_s18) 是唯一双受体（KDR+Flt-1）广谱阻断型**，且结合于二聚化核心——须以**双受体 SPR/BLI** 同时验证 KDR 与 Flt-1 阻断。
- **1✱ (#1_T0.3_s132) 计算结合最强**（ΔG 最负），但靶向 Flt-1+二聚面；须双受体 SPR 确认是否交叉阻断 KDR（否则体内效价存疑）。
- **全长 VEGF165 复核结论**：三主线候选在含 HBD 全长抗原下结合仍强且 vdW 全负（无穿模），设计稳健；相对排序一致（1✱ ≳ 4✱ > 2✱）。所有结合/功能实验**必须用全长 VEGF165**（见 §7）。

---

## 2. 合成 / 表达策略

**长度现实**：SPPS（固相多肽合成）经济上限约 50 aa。本批 60–100 aa，建议：

| 候选 | 长度 | 推荐方式 |
|---|---|---|
| #3 | 60 | SPPS（临界）或重组均可 |
| 1✱, 4✱ | 80 | **重组表达**（首选） |
| 2✱, #5 | 100 | **重组表达**（首选） |

**推荐重组路线**（E. coli BL21(DE3)）：
- N 端融合 **可切除标签** `His6–SUMO`（或 His6–Trx/MBP）提升可溶性；SUMO/TEV 蛋白酶切去除，得天然 N 端。
- 表达后：Ni-NTA → SEC-HPLC（Superdex 75）→ 单体峰收集。
- **二硫键处理**：最终构建体 Cys 数（以 `cro_order/constructs_cro.csv` 为准）：**1✱ #1_T0.3_s132：2 个（7,75）；2✱ #2c_T0.1_s18：2 个（18,86）；#3 vegf_len60_1：2 个（19,38）；4✱ T0.1_s35：0 个；NEG：0 个**。模型(Boltz/ESMFold)未强制建模二硫键，故结构为还原态几何；表达后经 **MS + 非还原 SDS-PAGE** 判定氧化态。多 Cys 候选用氧化折叠菌株；若 Cys 非界面关键，可 Cys→Ser/Ala 去聚集（见 `constructs/CLONING.md` 第 5 节）。
- **末端修饰**（可选）：N 端乙酰化 + C 端酰胺化提升蛋白酶稳定性。
- **浓度定量**：#1/#3 无 Trp，A280 不可用 → 用 AAA（氨基酸分析）或定量肽图；其余用 A280（Extinction Calculator 按序列计算）。

**对照**：
- 阴性：同长度同组成的 **scramble** 序列（或无关 de novo 微蛋白），验证特异性。
- 阳性：**bevacizumab (Avastin)**、**aflibercept (Eylea)** 作 VEGF165 结合/功能阳性对照。

---

## 3. 二级结构验证（圆二色谱 CD）

目的：确认候选**已折叠**（而非随机卷曲），并与预测结构一致。

- 缓冲液：低紫外吸收，如 10 mM NaPi pH 7.4 + 20 mM NaCl（避免 Tris/HEPES 紫外吸收）。
- 样品：0.1–0.5 mg/mL，光程 0.1 mm（或 1 mm 配稀释），**25 °C**，波长 **190–260 nm**，步长 1 nm。
- 热稳定性：222 nm 监测，20→90 °C 升温（2 °C/min），得 **Tm**。
- 判读：
  - 有特征性远紫外谱（~208/222 nm 双负峰提示 α-螺旋；~218 nm 单负峰提示 β-折叠）→ 已折叠。
  - 与 Boltz/AF3 预测二级结构含量比对（预测为折叠 miniprotein）。
  - 若呈随机卷曲 → 未折叠，需优化表达/复性或重新设计。
- 注：#1/#3 无 Trp，芳香族信号弱，但 CD 对主链敏感，不受影响。

---

## 4. 结合实验（SPR / BLI）— 核心验证

**中心目标**：(a) 确认 binder 结合 VEGF；(b) 量化 KD/kon/koff；(c) **判定是否阻断 VEGF–受体**（KDR 与 Flt-1 双靶点）。

### 4.1 直接结合（初筛，5 个候选）
- **BLI（Octet）优先初筛**：AHC/SSA 生物层固定 **VEGF165**（全长，关键点见 §7），浸入 binder 梯度（如 1 nM–5 µM，3 倍稀释，8–10 点），实时拟合 1:1 结合动力学。省样品、高通量。
- **SPR（Biacore）复核**：CM5 芯片 amine coupling 固定 VEGF165；binder 作为分析物梯度流过；多循环动力学（MCK）得 KD/kon/koff。固定 BSA 作背景对照。

### 4.2 双受体阻断实验（关键，对应结合模式结论）
在独立 SPR 通道分别固定 **KDR D2D3** 与 **Flt-1 D2**：
- **竞争法**：VEGF165 与 binder 预孵育（摩尔比梯度）后流经受体芯片；若结合信号被抑制 → binder 阻断该受体。
- 判定矩阵：

| 候选 | 预期阻断 KDR? | 预期阻断 Flt-1? | 验证意义 |
|---|---|---|---|
| 4✱ (T0.1_s35, KDR 位点) | 是（命中KDR位点） | 可能部分 | 主效价靶点确认 |
| 1✱ (#1_T0.3_s132, Flt-1+二聚面) | 待验证 | 是 | 须证明确实交叉阻断 KDR，否则效价存疑 |
| #3 (vegf_len60_1, Flt-1+二聚面) | 待验证 | 是 | Flt-1 机制探针 |
| 2✱ (#2c_T0.1_s18, 双受体) | 是（KDR+Flt-1） | 是 | **双受体广谱阻断**，关键验证 |
| #5 (vegf_len100_3, Flt-1主导) | 否（KDR 0） | 是 | 封闭 Flt-1；不阻断 KDR 信号，作 Flt-1 选择性探针 |
| NEG (scrambled) | 否 | 否 | 排除非特异性假阳 |

- 阳性对照：bevacizumab / aflibercept 对 VEGF165 的结合与竞争；阴性：scramble binder、无关微蛋白。
- 再生：10 mM Gly pH 2.0 或 2 M MgCl2，需优化避免 VEGF165 失活（勿用强酸长时）。

---

## 5. 功能 / 竞争细胞实验

- **HUVEC 增殖/存活**：VEGF165 刺激 + binder 梯度，CCK-8 / EdU，求 IC50。
- **管形成（tube formation）**：Matrigel 上 HUVEC + VEGF165 ± binder。
- **受体磷酸化 Western**：p-VEGFR2(Tyr1175) / p-VEGFR1(Tyr1333) —— **直接判定阻断 KDR 还是 Flt-1 信号**，与 §4.2 互为印证。
- 报告细胞：稳定表达 KDR 或 Flt-1 的工程细胞系（如有）。

---

## 6. 决策标准（Go / No-Go）

推进 lead optimization 的门槛（建议）：
- SPR/BLI **KD < 1 µM**（理想 < 100 nM）；
- 至少阻断 **KDR**（或 KDR+Flt-1 双阻断）；
- 细胞 **IC50 < 10 µM** 且 p-VEGFR 抑制明确；
- 单体、可溶、Tm > 40 °C。

按优先级：**4✱ 先行（KDR 位点 T0.1_s35，药效学最相关）→ 1✱/2✱ 并行（结合最强 + 双受体，须过双受体关）→ #3 跟进（Flt-1 探针）→ #5 备选（Flt-1 选择性）；全程含 NEG 阴性对照**。

---

## 7. 风险与注意

1. **模型仅含 VEGF 1–98（缺 C 端肝素结合区 HBD 111–165）**：binder 仅作用于 N 端受体结合域。**必须用全长 VEGF165**（而非本模型截短体）做所有结合/功能实验，否则可能漏掉 HBD 介导的 KDR D1 / αvβ3 信号未被阻断的风险（Takada 2024）。
2. **计算 ΔG 为相对富集值，非绝对亲和力**；最终以 SPR/BLI 的 KD 为准。
3. **穿模伪结构已系统剔除（vdW 双门控）**：早期个别候选（如旧 #2 de-A 版）初算 MM-GBSA 出现 vdW>0 穿模伪值（ΔG 数千 kcal），经 PDBFixer+OpenMM Cα 约束全原子最小化解除界面互穿后重算为真实负值。最终三主线候选（1✱/2✱/4✱）在**全长 VEGF165 复核**中 7/9 seed vdW 全负、另 2 个穿模 seed 已剔除（详见 `results/deA_redesign_1_2/report_fullVEGF_recheck.md`）。`structures/` 下 `.complex.pdb` 为 Boltz 预测结合复合物参考。
4. **可开发性未评分**：MW 均 <15 kDa 满足，但需补 T 细胞表位、聚集倾向（SEC/DLS）、蛋白酶稳定性评估。
5. **免疫原性/种属**：若用鼠源体系，注意人 VEGF-A 与鼠 Vegfa 同源度高（~88%），但受体交叉结合需实验确认。

---

## 8. 物料与订单速查

| 物料 | 用途 | 备注 |
|---|---|---|
| VEGF165（重组人，全长） | CD/SPR/BLI/功能 固定相或分析物 | 必须全长（含 HBD） |
| VEGF110（1–110 截短） | 对照（无 HBD） | 验证 HBD 依赖性 |
| KDR D2D3（胞外域） | 双受体阻断实验 | Sino Biological / 自制 |
| Flt-1 D2（胞外域） | 双受体阻断实验 | 同上 |
| bevacizumab / aflibercept | 阳性对照 | 市售 |
| scramble 微蛋白 | 阴性对照 | 同长度同组成重排 |
| Biacore / Octet 系统 | 动力学 | — |
| CD 光谱仪 (190–260 nm) | 二级结构 | 低紫外缓冲 |
