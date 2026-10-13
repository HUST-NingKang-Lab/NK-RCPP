# 图表与术语对齐说明 — NK-RCPP

> 权威来源：**NK-RCPP: a receptor-conditioned computational workflow for prioritizing neurokinin-receptor
> peptide candidates**（`paper_v5.docx`）。本文件只做对齐说明，凡与稿件冲突处一律以稿件为准。
> 图号 ↔ 文件 ↔ 生成脚本的完整对照见 [`../figures/README.md`](../figures/README.md)。

## 0. 论文定位（决定一切表述）

**回顾性、描述性审计**（retrospective, descriptive audit），不是设计成功故事，不是新的生成算法，不含湿实验。
稿件原文：

> *"This work was a retrospective, descriptive audit of a fixed computational screening repository… the
> name does not imply a new generation algorithm or experimental validation."*

> *"These findings quantify changes in selection and their uncertainty, rather than demonstrating improved
> experimental hit rates."*

受体重层审计**不计算 p 值**——阈值为确定性优先级规则，而非假设检验。探索性、未校正的 p 值只出现在外部
NK2R 基准（Fig. S6/S7）。

### 术语对照表（图内、正文必须使用右列）

| ❌ 禁用（早期提纲/旧图） | ✅ 稿件用语 |
|---|---|
| ΔPos | **Delta-on** = candidate score − receptor-matched control score |
| frac_safer、"off-target safety"、"脱靶安全性" | **fraction with a higher structural-neighbour score**（源字段仍叫 `frac_safer`） |
| off-target proteins、"脱靶蛋白"、"脱靶风险" | **structural neighbours**（*not validated pharmacological off-targets*） |
| ref2015 score | **REF15 total model score**（*not binding free energies*） |
| known ligand | **receptor-matched control** |
| dual-optimal / both / "双优" | **primary rule**（Delta-on < 0 且 fraction > 0.50 且 count proxy ≥ 0.90） |
| — | **strict sensitivity rule** ≥ 0.75（敏感性阈值，非独立验证的决策边界） |
| — | **legacy audit rule**（前两条判据，不含 count proxy） |
| — | **count proxy** = `n_common / Table S3 行数` |

> structural-neighbour fraction **不是** biological safety probability、不是 selectivity ratio；
> 凡 "safety"/"selectivity"/"安全性"/"选择性" 表述均超范围。

### 稿件声明的边界

1. REF15 total model score ≠ 结合自由能；受体系统之间**不可比**，只用受体内部 Delta-on。
2. structural neighbours 是异质集合（含受体结构本身、无关酶），**非**经验证药理脱靶。
3. 所审计的软件包中 pair-level 行缺失 → 聚合值原为**继承结果**；本仓库保留逐复合物打分，可重算
   （见 [`../step5_analysis/README.md`](../step5_analysis/README.md)）。
4. NK2R 14/32 候选缺聚合值 → primary 计数边界 10–20、strict 边界 2–12；该边界**不预测**缺失结果。
5. 阈值敏感：cutoff 由 > 0.50 提到 > 0.55 会移除**全部** NK1R 候选；阈值附近的候选互为
   *neighbouring priorities*，不是不同的生物学类别。

## 1. 关键数字速查（均与稿件一致）

```
候选总数 96（每受体 32）+ 3 个受体匹配对照
Delta-on < 0：NK1R 14/32（43.8%）| NK2R 23/32（71.9%）| NK3R 6/32（18.8%）
最优（score，Delta-on）：NK1R-18 517.29（−57.10）| NK2R-27 636.35（−134.78）| NK3R-21 1965.25（−275.29）
对照分数：Substance P 574.39 | NKA 衍生肽 771.13 | Neurokinin B 2240.54
近邻汇总覆盖：NK1R 32/32、NK2R 18/32、NK3R 32/32；合计 82 条可评估候选、4,518 次比较
面板（Table S3 行数）：66 / 27 / 60；所有可评估候选 count proxy 均为 1.00
primary rule：NK1R 3（NK1R-21、-29、-13）| NK2R 10（可评估）| NK3R 3（NK3R-21、-13、-19）
strict ≥0.75：仅 NK2R-04、NK2R-19（各 21/27 = 0.778；Delta-on −100.60、−72.21）
保留比例：NK1R 36/66、35/66、34/66 | NK3R 34/60、34/60、33/60
on-target 占优池 43 → 可评估 33 → 联合保留 16（缩减 51.5%）
等规模交集（3/10/3）：0/3、7/10、1/3 → 合计 8/16
不一致案例：NK1R-18 29/66 = 0.439 | NK2R-27 11/27 = 0.407 | NK1R-14 41/66 = 0.621（Delta-on +50.23）
共识序列：NK1R SVPRQQCFGF | NK2R PKPRQCCFQR | NK3R DTPYQQSNGF（NK2R-02 为 11-mer，位点分析剔除）
基准（n = 10）：z < 0 9/10 | Delta-on < 0 5/10 | fraction > 0.50 8/10 | 联合 3/10
  研究侧：18 条可评估中 14 条 fraction > 0.50、10 条联合；11/31 条 10-mer z < 0
  Spearman rho = +0.42（P = 0.23）；MWU P = 0.2204（原始分，10 vs 32）、0.1699（fraction，10 vs 18）、
  5.09 × 10⁻⁴（长度校正 z，10 vs 31）——全部探索性、双侧、未校正
序列空间：1,048 条（949 排行榜 + 99 研究）；92 条 7–10 残基构成小岛；4/10 竞赛肽在小岛内
  最近参考中位距离：10 条 primary NK2R 候选 0.04 vs 全部 99 条 0.07（二维投影单位）
```

## 2. 图号 ↔ 内容（文件见 `figures/`）

| 编号 | 内容要点 |
|:--|:--|
| **Fig. 1** | 流程 a–f：生成 → 靶内对接 → 结构近邻面板 → 配对比较 → 规则 → 受体特异的优先级 |
| **Fig. 2** | 受体分层审计：Delta-on 分布、Delta-on × fraction（count proxy 编码、缺失叉号）、计数类别、每受体最优 primary 候选 |
| **Fig. 3** | 候选级 fraction 排名（星号 = strict）+ 继承的配对均值差值 |
| **Fig. 4** | 位点氨基酸频率，共享 0–0.60 色标，共识阈值 0.30 |
| **Fig. 5** | a–d 归档口袋视图；e–h ESM-2/UMAP 序列空间、长度、最近参考距离、组成 |
| **S1** | 完整靶内 REF15 排名 |
| **S2** | 逐候选–对照中位配对差值 |
| **S3** | 面板注释组成 |
| **S4** | 代表性保留构象 |
| **S5** | 配对比较（散点、fraction、分布、示意 trace） |
| **S6** | 探索性 NK2R 竞赛对比（长度校正、分布、z vs EC50、fraction、trace 叠加） |
| **S7** | 描述性基准汇总（达标计数、记分卡、z 分布、按长度） |
| **S8** | 扩展对接诊断（对照 vs 最优候选、标准化分布、Delta-on 全景、trace） |
| **S9** | 扩展序列频率汇总，0–0.50 色标（与 Fig. 4 的 0–0.60 区分） |

## 3. 写作规范

1. 肽段首次出现用代号 + 序列（如 *NK1R-18*，`SPRRQLGFGL`），此后只用代号；图内一律代号。
2. 数字精度：分数 1–2 位小数，百分比 1 位。
3. 引用风格：Nature 编号制。
4. 凡涉及 structural-neighbour 的指标，**不得**表述为 "off-target safety"、"selectivity"、"safety"；
   凡涉及 docking 打分，**不得**表述为 binding free energy 或 affinity。
5. 前瞻性验证的表述：等规模、受体特异的 testing budget；strict 的 2 条 NK2R 候选作为预设亚组；
   先确认 identity/purity/solubility/aggregation/proteolytic stability，再做跨三受体的结合与功能
   （钙流或 IP1、cAMP、β-arrestin、内化）测定。

## 4. 参考文献

以稿件 References 1–29 为准（速激肽系统综述、NK1R/NK3R 结构与 NK1R/NK2R 结构–功能、
PepMLM/RAPiDock/REF15/Foldseek/CD-HIT/ESM-2/UMAP、生成与 binder 设计基准、GPCR 肽设计评估等）。
