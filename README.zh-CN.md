# NK-RCPP

**Neurokinin Receptor–Conditioned Peptide Prioritization（神经激肽受体条件化多肽优先级筛选）**
*面向 NK1R / NK2R / NK3R 的 96 条候选多肽的回顾性、受体配对审计*

[![Python](https://img.shields.io/badge/python-3.12-3776ab)](https://www.python.org/)
[![pandas](https://img.shields.io/badge/pandas-3.0.1-150458)](https://pandas.pydata.org/)
[![License: MIT](https://img.shields.io/badge/code-MIT-blue.svg)](LICENSE)
[![Data: CC BY 4.0](https://img.shields.io/badge/data%2Ffigures-CC%20BY%204.0-lightgrey.svg)](LICENSE)

![图 1 | NK-RCPP 工作流](figures/main/Fig1_workflow.png)

---

## 范围界定

本仓库是 NK-RCPP 的归档数据集、代码与图源。它被设计为对一份固定计算筛选仓库的**回顾性、描述性审计**：
只重新分析已归档的输出，**不**重跑多肽生成与对接。"NK-RCPP"指代被重建的流程（受体条件化生成 → 靶内对接 →
结构近邻聚合 → 确定性优先级规则），该名称不意味着新的生成算法或实验验证。

阅读时请牢记三条边界：

* **不含湿实验数据。** 这里没有任何证据表明被优先排名的多肽能够结合、激活或选择性作用于神经激肽受体。
  本审计量化的是：当配对证据、缺失情况与判定阈值被显式写出时，候选清单**如何改变**；它并不宣称更高的
  命中率、亲和力、激动活性、选择性或安全性。
* **REF15 总模型打分不是结合自由能。** 它依赖所建模的复合物，因此比较只限于同一受体体系内的候选与
  受体匹配对照之间。
* **结构近邻不是经实验验证的药理学脱靶，** 近邻打分比例不是安全概率、也不是选择性比值；它是对打分行为
  在配对情境下的一次计算压力测试。

受体层面的审计**未计算任何 p 值**——其阈值是确定性优先级规则，而非假设检验。探索性、未校正的 p 值仅出现在
第 6 步的外部 NK2R 基准中。

## 主要结果

| | NK1R | NK2R | NK3R |
|:--|--:|--:|--:|
| 生成的候选数（PepMLM） | 32 | 32 | 32 |
| 受体匹配对照 | P 物质 | NKA 衍生肽 | 神经激肽 B |
| 对照 REF15 总模型打分 | 574.39 | 771.13 | 2240.54 |
| Delta-on < 0 的候选 | **14/32**（43.8%） | **23/32**（71.9%） | **6/32**（18.8%） |
| 最优候选（打分，Delta-on） | NK1R-18（517.29，−57.10） | NK2R-27（636.35，−134.78） | NK3R-21（1965.25，−275.29） |
| 具有近邻汇总的候选 | 32/32 | 18/32 | 32/32 |
| 结构近邻面板（Table S3 行数） | 66 | 27 | 60 |
| 主规则候选数 | 3 | 10 | 3 |
| 严格规则（≥ 0.75）候选数 | 0 | 2 | 0 |

三个受体合计：**82 条可评估候选**、**4,518 次配对结构近邻比较**。在 33 条可评估的"靶内打分占优"候选中，
联合主规则保留了 **16 条**，将该候选池**缩减 51.5%**。在同等候选清单规模下，联合选择与"仅靶内打分"选择
仅共享 **8/16** 条——即使测试预算不变，**一半的入选身份发生了改变**。

> **主规则** —— Delta-on < 0 **且** 近邻打分比例 > 0.50 **且** 报告比较数覆盖代理 ≥ 0.90，在各受体内部独立施加。
> **严格敏感性规则** —— 同上，但比例阈值取 ≥ 0.75。
> **旧版审计规则** —— 仅前两条，不含覆盖代理。

## 流程

| 步骤 | 操作 | 工具 | 目录 |
|:--|:--|:--|:--|
| 1 | 受体条件化候选生成：每受体 32 条多肽（多为 10 残基）+ 受体匹配对照 | PepMLM（作者公开 Colab） | [`step1_peptide_generation/`](step1_peptide_generation/) |
| 2 | 靶内对接：每条多肽 10 个构象，保留最低 REF15 总模型打分；Delta-on = 候选 − 受体匹配对照 | RAPiDock + Rosetta REF15 | [`step2_positive_docking/`](step2_positive_docking/) |
| 3 | 结构近邻检索：PDB100 的 3Di 与氨基酸联合检索，随后去冗余 | Foldseek、CD-HIT 0.50 | [`step3_foldseek_search/`](step3_foldseek_search/) |
| 4 | 配对结构近邻对接：候选**与**对照都对接至每个面板成员，保证比较是配对的 | RAPiDock + REF15 | [`step4_negative_docking/`](step4_negative_docking/) |
| 5 | 聚合与优先级：近邻比例、覆盖代理、主/严格规则、阈值网格、缺失边界、清单交集 | Python（pandas/NumPy） | [`step5_analysis/`](step5_analysis/) |
| 6 | 探索性外部 NK2R 基准：10 条竞赛肽、长度匹配随机对照、序列空间投影 | Python；归档的 ESM-2 + UMAP 坐标 | [`step6_competition_benchmark/`](step6_competition_benchmark/) |

图片采用论文编号，[`figures/README.md`](figures/README.md) 给出每个文件对应的 `Fig.`/`Supplementary Fig.`
编号与生成脚本。

## 目录结构

```
NK-RCPP/
├── step1_peptide_generation/   input/  受体链序列            output/ 候选表 + 伪困惑度
├── step2_positive_docking/     input/  受体 PDB、RAPiDock 输入 CSV   output/ REF15 打分（长表 + 宽表）
├── step3_foldseek_search/      output/ 各受体 PDB100 命中、比对、面板 FASTA
├── step4_negative_docking/     input/  按批次切分的对接输入   output/ 逐复合物近邻打分
├── step5_analysis/output/      排名、配对聚合、逐对差值、UMAP 坐标
├── step6_competition_benchmark/output/ 竞赛肽与长度对照打分
├── figures/                    main/ supplementary/ tables/（Fig. 1–5、S1–S9、TableS1–S4）
├── scripts/                    生成全部表格与图片的脚本
└── docs/                       论文指南、竞赛对比、项目状态
```

## 归档中的结果

### 靶内打分

靶内对接打分优于匹配对照的候选比例以 NK2R 最高（23/32，71.9%），其次 NK1R（14/32，43.8%）与
NK3R（6/32，18.8%）。各受体内部改善幅度最大者为 NK1R-18（517.29 对 574.39，Delta-on = −57.10）、
NK2R-27（636.35 对 771.13，Delta-on = −134.78）、NK3R-21（1965.25 对 2240.54，Delta-on = −275.29）。
NK1R 的比较还将 10 残基候选与 11 残基的 P 物质配对。

### 联合主规则

匹配的结构近邻汇总覆盖 **82** 条候选、**4,518** 次比较；面板分别含 NK1R 66、NK2R 27、NK3R 60 条注释行，
且所有可评估候选的报告覆盖代理均为 **1.00**，故 ≥ 0.90 的代理条件未进一步区分候选。

| 受体 | 主规则清单 | 近邻打分比例 | Delta-on |
|:--|:--|:--|--:|
| NK1R | NK1R-21、NK1R-29、NK1R-13 | 36/66、35/66、34/66 | −11.57、−29.80、−15.70 |
| NK2R | NK2R-04、NK2R-19、NK2R-14、NK2R-07、NK2R-26、NK2R-09、NK2R-11、NK2R-10、NK2R-13、NK2R-22 | NK2R-04 与 NK2R-19：21/27（0.778） | −100.60、−72.21 |
| NK3R | NK3R-21、NK3R-13、NK3R-19 | 34/60、34/60、33/60 | −275.29、−55.26、−28.47 |

仅 **NK2R-04** 与 **NK2R-19** 同时满足严格规则（≥ 0.75）。

### 两个排序维度之间的不一致

靶内排序与结构近邻排序反映打分行为的不同侧面，最优候选并不重合：

| 候选 | 靶内排名 | 近邻打分比例 | 结果 |
|:--|:--|:--|:--|
| NK1R-18 | NK1R 最优（Delta-on −57.10） | 29/66 = 0.439 | 未通过比例条件 |
| NK2R-27 | NK2R 最优（Delta-on −134.78） | 11/27 = 0.407 | 未通过比例条件 |
| NK1R-14 | Delta-on **+50.23** | 41/66 = 0.621（合格 NK1R 中最高） | 未通过靶内条件 |

两个终点承载不同信息，因此保持彼此独立；联合使用可降低对单一打分的依赖。

### 阈值与缺失

* 将比例阈值由 > 0.50 提高到 **> 0.55** 会**移除全部 NK1R 候选**。
* 严格规则 **≥ 0.75** 仅保留上述两条 NK2R 候选。
* 这些变化源于边界附近的离散计数——阈值上下相邻的候选应视为**相邻优先级，而非不同的生物学类别**。
* **14 条 NK2R 候选**缺少结构近邻汇总，其中 10 条 Delta-on < 0。因此无约束边界为：主规则
  **10–20 条**、严格规则 **2–12 条**。该边界量化的是计数对不可得证据的依赖，并不预测缺失结果。

### 序列模式

在每位点频率阈值 0.30 下，共识序列为 **SVPRQQCFGF**（NK1R）、**PKPRQCCFQR**（NK2R）、
**DTPYQQSNGF**（NK3R）。NK2R-02 含 11 个残基，从位点热图中剔除（NK2R 余 31 条），但仍保留在候选清单与
原始靶内分析中。这些分布描述的是固定的生成集合，可用于多样性感知的实验抽样；在缺少富集与重复生成证据时，
**不能称为受体结合基序**。

### 探索性 NK2R 基准（第 6 步）

| 终点 | 竞赛肽（n = 10） | 研究中的 NK2R 候选 |
|:--|:--|:--|
| 长度校正 z < 0 | 9/10 | 11/31 条 10 残基候选 |
| Delta-on < 0 | 5/10 | — |
| 近邻打分比例 > 0.50 | 8/10 | 14/18 条可评估候选 |
| 联合规则 | 3/10（竞赛 9 肽中 3/5） | 10/18 条可评估候选 |

打分与效价的关联不确定（Spearman rho = +0.42，P = 0.23，n = 10）。探索性双侧 Mann-Whitney 检验：
原始打分为 P = 0.2204（10 对 32 条肽）、近邻比例为 P = 0.1699（10 对 18）、长度校正 z 为
P = 5.09 × 10⁻⁴（10 对 31）。这些比较**未做多重性校正**，且涉及不同长度、不同组成基线与非同一抽样框，
不能据此宣称结合或活性更优。

### 序列空间

归档的 ESM-2 + UMAP 投影共 **1,048** 条肽：949 条仅来自排行榜，99 条来自本研究（96 条候选 + 3 条对照）。
其中 92 条 7–10 残基序列构成分离的小岛，956 条 3–40 残基序列位于主投影；10 条竞赛肽中有 4 条落在小岛内。
到最近保留排行榜条目的中位距离：10 条主规则 NK2R 候选为 **0.04**，全部 99 条研究条目为 **0.07**。
距离以任意二维投影为单位，**不是**原始 ESM-2 表示空间，也不构成新颖性、功能或受体特异性的证据。

## 术语对照：论文 ↔ 归档

论文对若干归档字段作了更精确的改名，阅读 CSV 时请对照：

| 论文术语 | 归档字段 | 定义 |
|:--|:--|:--|
| Delta-on | *派生量* | 候选总打分 − 受体匹配对照总打分；**< 0 更优** |
| REF15 总模型打分 | `best_score` / `rank1_score` | 该建模复合物可得的 `ref2015score` 最小值；不是结合自由能 |
| 结构近邻打分更高的比例 | `frac_safer` | `n_safer_delta_pos / n_common` |
| 打分更高的近邻数 | `n_safer_delta_pos` | 候选打分高于对照打分的比较次数 |
| 共同比较数 | `n_common` | 候选与对照均被打分的结构近邻数 |
| 报告比较数覆盖代理 | *派生量* | `n_common / 该受体 Table S3 行数` |
| 结构近邻 | TableS3 行 | Foldseek 检索得到的比较结构，**不是**经验证的脱靶 |

## 数据约定与坑

* **第 4 步 NK2R 的批次标签不是规范候选编号。** `nk2r_neg_inputs/peptide01…18` 对应
  NK2R-01/03/04/07/09/10/11/13/14/19/22/23/24/26/27/28/30/32 这一子集，而 `peptide19` 是**对照**。
  NK1R 与 NK3R 则对齐（`peptideNN` 即 `NK1R-NN` / `NK3R-NN`，对照为 `known`）。请一律按**序列**合并，
  绝不按批次标签——`scripts/paired_neighbour_deltas.py` 正是这样做的。
* **旧版审计规则与主规则的区别**只在于覆盖代理；在本归档中两者结果一致，因为可评估候选的代理均为 1.00。
* **`rank1` 与 `best_score`。** 下游分析一律使用 `rank1`（RAPiDock 排名第一的构象）；
  `best_score`（10 个构象中的最小值）一并提供以便对比。
* **`COMP_TOP10_negative_docking_input.STALE-offby-one.csv` 是刻意保留的。** 该合并版竞赛输入文件的标签
  相对逐批次输入整体错位一格，按复合物名查序列会静默取到相邻多肽；以 `comp_top10_neg_inputs/` 下的
  逐批次文件为准，`scripts/harvest_docking_scores.py` 会在不一致时告警。
* **NK2R 伪困惑度为空。** NK2R 的 Colab 导出未归档，序列由对接输入文件恢复；所有受体都无法重建原生生成器排名。
* **任何地方都不跨受体比较原始 REF15 打分**——其绝对值依赖所建模的体系。

## 关于证据层级

论文所审计的软件包只有聚合值：逐对结构近邻打分行缺失，因此候选级聚合只能作为"继承结果"处理。
**本归档比那个包更完整**——第 4 步提供底层的逐复合物 REF15 打分。运行

```bash
python scripts/paired_neighbour_deltas.py --repo .
```

即可从原始数字重建每个聚合背后的逐对比较，并精确复现随包发布的计数：

```
[ok] 4518 pair-level comparisons over 82 candidates
[check] NK1R: 32 candidates compared, 0 disagree with the shipped aggregates
[check] NK2R: 18 candidates compared, 0 disagree with the shipped aggregates
[check] NK3R: 32 candidates compared, 0 disagree with the shipped aggregates
```

因此，论文引用的 4,518 次比较 / 82 条候选在本仓库中是可复算的，而不只是继承值。仍然无法重建的是原始归档中
本就缺失的部分：RAPiDock 原始构象目录、原始生成输出，以及 COMP/长度对照作战的逐对行。

## 复现

```bash
# 审计算术（Python 3.12；论文使用 pandas 3.0.1 / NumPy 2.3.5）
python scripts/paired_neighbour_deltas.py --repo .     # 82 条候选、4,518 次逐对比较
python scripts/make_candidate_tables.py --repo .       # 第 1 步候选表
python scripts/harvest_docking_scores.py --repo .      # 由原始运行目录汇总第 2/4/6 步打分表
python scripts/harvest_foldseek_hits.py --repo .       # 第 3 步近邻面板
python scripts/make_peptide_codes.py --repo .          # TableS1
```

以上脚本**都不重跑**对接；对接所用 RAPiDock 配置见
[`step2_positive_docking/README.md`](step2_positive_docking/README.md)（`N = 10`、`inference_steps` /
`actual_steps = 16`、`ref2015`、本地权重 `rapidock_local.pt`、`conformation_partial = 1:1:1`、
`batch_size = 4`），GPU 运行脚本位于 `scripts/`。

图片重绘清单见 [`figures/README.md`](figures/README.md)。

## 局限

单一归档数据集；回顾性、未校准的打分阈值；无湿实验验证。96 条候选的原生 PepMLM 排名与活性数据不可得，
因此清单压缩无法证明命中率提升或结合物保留。缺失的生成设置、原始构象、逐对近邻打分与化学制备记录，
限制了对末端修饰、二硫键与面板组成的复现与评估。10 肽基准及每个长度 5 个组成特异随机对照只提供情境，
而非独立的药理学验证。独立的阈值校准、正交结构评估，以及规模相等的前瞻性结合与功能比较，是下一步的检验。

## 引用

若使用本归档，请引用 NK-RCPP 论文（正式发表后补全引用信息）及其所依赖的工具：
PepMLM（[Chen et al., *Nat. Biotechnol.* 2026](https://doi.org/10.1038/s41587-025-02761-2)）、
RAPiDock（[Zhao et al., *Nat. Mach. Intell.* 2025](https://doi.org/10.1038/s42256-025-01077-9)）、
Rosetta REF15（[Alford et al., *JCTC* 2017](https://doi.org/10.1021/acs.jctc.7b00125)）、
Foldseek（[van Kempen et al., *Nat. Biotechnol.* 2024](https://doi.org/10.1038/s41587-023-01773-0)）、
CD-HIT（[Li & Godzik, *Bioinformatics* 2006](https://doi.org/10.1093/bioinformatics/btl158)）、
ESM-2（[Lin et al., *Science* 2023](https://doi.org/10.1126/science.ade2574)）与
UMAP（[McInnes et al., *JOSS* 2018](https://doi.org/10.21105/joss.00861)）。

## 许可

代码采用 MIT 许可；数据、表格与图片采用 CC BY 4.0，详见 [`LICENSE`](LICENSE)。
NK2R 受体模型源自 AlphaFold2 预测，见 [`docs/ALPHAFOLD_TERMS_OF_USE.md`](docs/ALPHAFOLD_TERMS_OF_USE.md)。

---

<sub>English documentation: [`README.md`](README.md).</sub>
