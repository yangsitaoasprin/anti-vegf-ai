# #1–#5 重组表达构建体与引物设计

## 1. N 端融合策略 (实验设计 §2 约定)
- 融合：`His6 – Smt3(SUMO, 酵母) – [SUMO蛋白酶切] – binder(天然N端)`
- 表达后 **SUMO 蛋白酶 (Ulp1/SENP 家族)** 在 Smt3 C 端 Gly 后切断，释放**天然 N 端 binder**，无多余残基。
- 备选：**His6 – TEVsite(ENLYFQG) – binder**，用 TEV 蛋白酶切（若偏好 TEV 体系）。
- 纯化：Ni-NTA（His6）→ SEC-HPLC（Superdex 75）单体峰。

## 2. 密码子优化
- `constructs_dna.fasta` 为 E. coli 高表达偏向密码子序列（模板）。
- **下单前建议用专业工具二次优化**（GenSmart / IDT / Eurofins），并去除重复序列(SSR)、弱化 Shine-Dalgarno 样序列与发夹。

## 3. 克隆方案
### 3.1 Gibson 组装（推荐，无内部酶切位点限制）
- 将 binder 合成片段（含两端 20–24 nt 与线性化 pET-SUMO 载体同源臂）组装入载体。
- 载体同源臂需按**实际所用 pET-SUMO 载体图谱**替换（本表 `GIBSON_FWD/REV` 中的 18 nt 臂为示意，须替换为真实载体序列）。
### 3.2 限制酶克隆 (NdeI / XhoI)
- Fwd 引物含 NdeI(CATATG)+起始 ATG，Rev 含 XhoI(CTCGAG)+终止密码子（与 CRO 投递单一致）。
- 若 binder 内部出现 NdeI/XhoI 位点（见 `note` 列），须做**沉默突变**或改用 Gibson。

## 4. 引物表（与 constructs_cro.csv 完全一致，NdeI/XhoI + Gibson）

| 候选 | 构造ID | NdeI/XhoI FWD | NdeI/XhoI REV | Gibson FWD | Gibson REV | 备注 |
|---|---|---|---|---|---|---|
| #1 | #1_T0.3_s132 | `AAGAAGGAGATATACATATGCATCATCATCATCATCATAG` | `CCGCTCGAGTTACGCTTCCATCGCCTGGCATTC` | `GGGAATTCCATATGCATCATCATCATCATCATAGC` | `CGCTTCCATCGCCTGGCATTTGGTGGTGGTGGTGCTCGAG` | 无内部NdeI/XhoI位点, 两种克隆均可 |
| #2 | #2c_T0.1_s18 | `AAGAAGGAGATATACATATGCATCATCATCATCATCATAG` | `CCGCTCGAGTTACAGCGCCGCCGCCTGCGCTTT` | `GGGAATTCCATATGCATCATCATCATCATCATAGC` | `CAGCGCCGCCGCCTGCGCTTTGGTGGTGGTGGTGCTCGAG` | 无内部NdeI/XhoI位点, 两种克隆均可 |
| #3 | vegf_len60_1_T0.2_s1 | `AAGAAGGAGATATACATATGCATCATCATCATCATCATAG` | `CCGCTCGAGTTAACGCAGCTGACGCACCAGTTC` | `GGGAATTCCATATGCATCATCATCATCATCATAGC` | `ACGCAGCTGACGCACCAGTTTGGTGGTGGTGGTGCTCGAG` | 无内部NdeI/XhoI位点, 两种克隆均可 |
| #4 | T0.1_s35 | `AAGAAGGAGATATACATATGCATCATCATCATCATCATAG` | `CCGCTCGAGTTAACGTTTCGCCGCCACCGCCGC` | `GGGAATTCCATATGCATCATCATCATCATCATAGC` | `ACGTTTCGCCGCCACCGCCGTGGTGGTGGTGGTGCTCGAG` | 无内部NdeI/XhoI位点, 两种克隆均可 |
| #5 | vegf_len100_3_T0.1_s2 | `AAGAAGGAGATATACATATGCATCATCATCATCATCATAG` | `CCGCTCGAGTTATTTCGCCAGTTTCGCCGCTTC` | `GGGAATTCCATATGCATCATCATCATCATCATAGC` | `TTTCGCCAGTTTCGCCGCTTTGGTGGTGGTGGTGCTCGAG` | 无内部NdeI/XhoI位点, 两种克隆均可 |


## 5. 二硫键注意（Cys 数随重设计改变，以下为最终投递构造）

| 候选 | 构造ID | Cys 数 | 位置 | 提示 |
|---|---|---|---|---|
| #1 | #1_T0.3_s132 | 2 | 7,75 | 2 对潜在分子内二硫键（7,75）；建议氧化折叠菌株 (Origami / SHuffle / 周质) |
| #2 | #2c_T0.1_s18 | 2 | 18,86 | 2 Cys(18,86)；关注聚集 |
| #3 | vegf_len60_1_T0.2_s1 | 2 | 19,38 | 大概率 1 对二硫键（19,38） |
| #4 | T0.1_s35 | 0 | — | 0 Cys（去淀粉样替代版，无游离硫醇风险） |
| #5 | vegf_len100_3_T0.1_s2 | 1 | 92 | 1 Cys(92)，游离硫醇风险→聚集；如非界面关键可 Cys→Ser |

- 模型(Boltz/ESMFold)均**未强制建模二硫键**，结构中的 Cys 处于还原态几何；真实表达后须经 **MS + 非还原 SDS-PAGE + CD** 判定氧化态。
- 多 Cys(#1) 建议氧化折叠环境；**单 Cys 候选若不在结合界面且不参与结合，可评估 Cys→Ser/Ala 突变**，消除游离硫醇导致的聚集/二硫错配。
- 阴性对照 NEG_scrambled_T0.1_s35 为 0 Cys，无需二硫键处理。
