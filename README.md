# NK 受体激动剂短肽设计 — 结果与代码包

> **给接手 agent 的 README**：本包包含 NK1R/NK2R/NK3R 选择性激动剂短肽设计项目的**最终分析结果、图表与全部代码**，以及与 2025 多肽设计竞赛 Top10 作品的基准对比。不包含 RAPiDock 模型权重、原始对接构象等大文件。

## 1. 项目一句话

用 **PepMLM（条件肽生成）→ RAPiDock（扩散对接，ref2015 打分）→ Foldseek 3diaa（结构相似蛋白搜索）→ 负向对接（脱靶筛选）** 的端到端管线，为速激肽受体 NK1R/NK2R/NK3R 设计 10-mer 短肽激动剂候选；多个候选在**结合强度**与**脱靶选择性**上同时优于天然配体。并将同一流程应用于 2025 多肽设计竞赛排行榜前 10 名（有实测 EC50 的作品）做基准对比。

## 2. 目录结构

```
package/
├── README.md                          # 本文件
├── MANUSCRIPT_GUIDE.md               # 论文撰写指南（结构、数据引用、图表清单）
├── COMPETITION_COMPARISON.md         # 竞赛 Top10 对比分析（含长度校正方法学）
├── PROJECT_STATUS.md                 # 项目状态与路径/环境备忘
├── analysis/                         # 所有汇总数据 (CSV)
│   ├── {NK1R,NK2R,NK3R}_summary.csv              # 正/负向汇总
│   ├── {NK1R,NK2R,NK3R}_final_ranking.csv        # 最终排名 (ΔPos/ΔNeg/frac_safer)
│   ├── {NK1R,NK2R,NK3R}_pairwise_selectivity.csv # 逐蛋白配对 delta
│   ├── peptide_codes.csv                         # 代号↔序列映射
│   ├── competition_top10.csv                     # 竞赛 Top10 元数据 (含 EC50)
│   ├── COMP_TOP10_summary.csv                    # 竞赛肽对接结果
│   ├── NK2R_comparison.csv                       # 合并对比表
│   ├── length_control.csv                        # 长度对照实验原始数据
│   └── all_receptors_summary.csv
├── figures/                           # 论文图表 (SVG/PDF 矢量, 文本可编辑)
│   ├── Fig2_docking_performance.*     # 主图2: 正向对接 (4 panels)
│   ├── Fig3_selectivity.*             # 主图3: 脱靶选择性 (6 panels)
│   ├── Fig4_sequence_features.*       # 主图4: 序列特征 (4 panels)
│   ├── Fig5_competition_comparison.*  # 主图5: 竞赛对比细节 (6 panels)
│   ├── Fig6_head_to_head.*            # 主图6: 头对头总览 (4 panels: 达标率/计分卡/z分布/长度分组)
│   ├── Fig8_latent_space.*            # 主图8: ESM-2 潜在空间 UMAP + 4 个结合姿态 (8 panels)
│   ├── FigS1_positive_ranking.*       # 全候选正向排名
│   ├── FigS2_delta_neg_distribution.* # ΔNeg 排名
│   ├── FigS3_offtarget_composition.*  # 脱靶蛋白类型
│   ├── FigS4_binding_modes.*          # 结合模式 (6 panels)
│   └── tables/
│       ├── TableS1_peptide_codes.csv      # 代号→序列 (99条)
│       ├── TableS2_docking_data.csv       # 完整对接数据
│       ├── TableS3_offtarget_proteins.csv # 脱靶蛋白身份 (153条)
│       └── TableS4_competition_comparison.csv # 竞赛对比 (43条)
└── scripts/                           # 全部代码 (Python; 依赖 pandas numpy matplotlib scipy biopython)
    ├── summarize_results.py           # 1. 汇总正/负向分数 (按序列关联!)
    ├── pairwise_selectivity.py        # 2. 逐蛋白配对选择性
    ├── final_ranking.py               # 3. 最终排名
    ├── make_peptide_codes.py          # 4. 代号映射
    ├── make_fig2.py / make_fig3.py / make_fig4.py   # 5-7. 主图
    ├── make_fig5_comparison.py        # 8. 竞赛对比细节图
    ├── make_fig6_headtohead.py        # 9. 头对头总览图
    ├── make_fig8.py                   # 9b. ESM-2 潜在空间图 (UMAP + 结合姿态)
    ├── make_fig8_umap.py              # 9c. UMAP 降维
    ├── make_umap_data.py              # 9d. UMAP 输入数据准备
    ├── embed_esm.py                   # 9e. ESM2-650M 编码 (gpu3, RAPiDock env)
    ├── pdb_render.py                  # 9f. PDB 结合姿态渲染器 (纯 Bio.PDB+matplotlib)
    ├── make_supp_figs.py              # 10. 附图
    ├── make_binding_modes.py          # 11. 结构图
    ├── make_supp_tables.py            # 12. 附表
    ├── make_comp_top10_inputs.py      # 13. 竞赛 Top10 → 对接输入
    ├── make_length_control.py         # 14. 长度对照实验设计
    ├── summarize_comp_top10.py        # 15. 竞赛结果汇总 (含长度校正)
    ├── make_comp_table.py             # 16. 竞赛对比表 S4
    ├── foldseek_search.py             # Foldseek API (需联网)
    ├── parse_and_download.py          # 解析+CD-HIT(或Python聚类)+下载PDB
    ├── extract_chains.py              # 子链提取
    ├── cluster_sequences.py           # Python 聚类器 (cd-hit 替代)
    ├── create_neg_input.py            # 负向输入 CSV 生成
    ├── download_pdbs_parallel.sh      # 并行下载
    ├── run_docking_gpu3.sh            # GPU 对接命令 (gpu3 节点)
    ├── comp_top10_worker.sh           # 竞赛对接 worker (指定 GPU + 批次)
    └── fix_shifted_labels.py          # 标签错位修正脚本 (见 §6)
```

## 3. 复现分析流程

```bash
conda install pandas numpy matplotlib seaborn scipy biopython

# 数据汇总 (输入: peptide_docking/RAPiDock/results/ 下原始对接结果)
python scripts/summarize_results.py
python scripts/pairwise_selectivity.py
python scripts/final_ranking.py
python scripts/make_peptide_codes.py

# 图表生成 (输出: figures/)
python scripts/make_fig2.py && python scripts/make_fig3.py && python scripts/make_fig4.py
python scripts/make_fig5_comparison.py && python scripts/make_fig6_headtohead.py
python scripts/make_supp_figs.py && python scripts/make_binding_modes.py
python scripts/make_supp_tables.py

# Fig8 (ESM-2 潜在空间)
python scripts/make_umap_data.py            # 1. 组装 1048 条序列
# (在 gpu3 的 RAPiDock env 运行编码, ~2 分钟; 模型走本地缓存, 离线)
python scripts/embed_esm.py --input analysis/umap_input.csv --output analysis/esm_emb.npy
python scripts/make_fig8_umap.py            # 2. L2归一化 + UMAP 降维
python scripts/make_fig8.py                 # 3. 出图

# 竞赛对比
python scripts/make_comp_top10_inputs.py    # 生成对接输入
python scripts/make_length_control.py       # 长度对照输入
# (在 gpu3 运行对接: comp_top10_worker.sh)
python scripts/summarize_comp_top10.py
python scripts/make_comp_table.py
```

## 4. 关键数据口径（重要！）

- **ref2015 分数越低越好**（Rosetta 全原子能量，非结合自由能）
- **不同受体系统间绝对值不可比**，一律用受体内部 Δ（候选 − 已知）
- **跨肽长不可直接比较**：ref2015 约 +100.6 分/残基（R²=0.57）→ 必须用同长随机对照做 z 值校正
- **负向数据必须按序列关联**：NK2R 的 peptideNN 输入编号 ≠ 候选编号，按序列匹配才正确（脚本已处理）
- 指标定义：
  - `ΔPos` = 候选正分 − 已知正分（<0 = 结合目标更强）
  - `ΔNeg` = 候选负向均值 − 已知负向均值（>0 = 脱靶结合更弱）
  - `frac_safer` = 逐蛋白配对 delta>0 的脱靶蛋白比例（>0.5 更安全）
  - `both` = ΔPos<0 且 frac_safer>0.5（推荐筛选标准）
  - `z_vs_length` = (分数 − 同长随机对照均值) / 同长对照 SD
- **已知肽代号**：NK1R-K (RPKPQQFFGLM, Substance P)、NK2R-K (HKTDSFVGLM)、NK3R-K (DMHDFFVGLM, Neurokinin B)
- 图中肽段全部用代号（NK1R-01 等）；竞赛肽用 COMP-01..COMP-10，序列见 TableS1/S4

## 5. 核心结果摘要

**三个受体（本研究）**

| 受体 | 已知肽正分 | 双优候选 (ΔPos<0 且 frac_safer>0.5) |
|------|-----------|-------------------------------------|
| NK1R | 574.39 | NK1R-29, NK1R-13, NK1R-21 (3个) |
| NK2R | 771.13 | NK2R-04/07/09/10/11/13/14/19/22/26 (10个) |
| NK3R | 2240.54 | NK3R-21, NK3R-13, NK3R-19 (3个) |

- 共识基序：NK1R SVPRQQCFGF | NK2R PKPRQCCFQR | NK3R DTPYQQSNGF
- 脱靶蛋白：NK1R 66 个、NK2R 27 个、NK3R 60 个（多为 GPCR：H1、M2/3/5、阿片 μ/δ/κ、CXCR4、5-HT、CB1、PAR1/2 等）

**竞赛 Top10 基准对比（Fig5 + Fig6）**

| 指标 | 竞赛 Top10 | 本研究 (NK2R) |
|------|-----------|--------------|
| ΔPos < 0 | 5/10 (50%) | 23/32 (72%) |
| z_vs_length < 0 | 9/10 (90%) | 11/32 (34%) |
| frac_safer > 0.5 | 8/10 (80%) | 14/18 (78%) |
| **双优 (both)** | **3/10 (30%)** | **10/18 (56%)** |
| 同长度段双优 | 9-mer: 3/5 (60%) | 10-mer: 10/18 (56%) |

- 长度对照（随机组成匹配肽，n=5/长度）：6-mer 625.7±14.4 | 9-mer 780.0±105.9 | 10-mer 687.4±52.0 | 15-mer 1502.0±547.8
- 对接打分对实验 EC50 的预测力：ρ=+0.42 (P=0.23, n=10)，方向正确但无显著性
- 竞赛前 10 是高质量设计：9/10 显著优于同长随机对照（非长度假象）

## 6. 已知局限与数据陷阱（必读）

1. **NK2R 有 14/32 候选缺负向数据**（原始 valid_pep 子集选择），含正向最优 NK2R-21/16/08；TableS2 中 N_offtarget=0 即缺失
2. **负向对接标签曾错位（已修正）**：初次生成负向输入时编号 off-by-one（comp00–comp09 vs comp01–comp10），导致 1 号肽漏跑、其余肽标签前移一位。已通过 `fix_shifted_labels.py` 逆序重命名修正并补跑 comp01，最终 297/297 完整。**如重新生成输入，务必校验 compNN ↔ 序列的对应关系**
3. **亚型选择性在数据中不显著**：NK2R 候选对 NK1R 结构 (6j20) 的 ΔNeg ≈ +1.5~+17.9（NK2R-27 为 −8.5），与已知肽相当
4. **长度校对照样本量小**（每长度 n=5），15-mer 组 SD 极大（548），其 z 值被压缩；跨长度组比较 z 值需谨慎
5. 纯计算研究，无 wet-lab 验证；ref2015 非结合自由能；扩散对接受初始构象影响
6. Foldseek 3diaa 命中含糖基转移酶/合酶类（局部拓扑相似），非药理学相关脱靶
7. **6-mer CFYFLM (COMP-03) 是长度陷阱案例**：原始分数最低但 frac_safer 仅 0.19，说明单用原始分数会误判

## 6b. 图表样式（figures4papers house style）

所有图遵循 [figures4papers](https://github.com/ChenLiu-1996/figures4papers) 的发表级绘图规范，
统一定义在 `scripts/fig_style.py`：

```python
from fig_style import PALETTE, RECEPTOR_COLORS, apply_style, save_pub
apply_style()                       # 设定 house style rcParams
save_pub(fig, "figures/FigX")       # 导出 SVG + PDF + 600 dpi LZW-TIFF
```

采用的规范要点：
- **语义配色**：`PALETTE` 取自 `references/api.md`（blue_main `#0F4D92`、
  green_3 `#8BCF8B`、red_strong `#B64342`、teal `#42949E`、violet `#9A4D8E`、neutral `#CFCECE`）
- **受体配色**：NK1R=blue_main、NK2R=teal、NK3R=violet（`RECEPTOR_COLORS`）
- **加粗坐标轴与刻度**：`axes.linewidth = 1.1`（house style 比默认更粗）
- **无网格、无图例边框、仅左/下轴脊**
- **导出**：SVG（文本可编辑）+ PDF（TrueType 嵌入）+ 600 dpi LZW-TIFF

> ⚠️ 标定说明：figures4papers 的 demo 用大画布（28 in 宽 / font 24 / lw 3），
> 本项目按稿件图幅 183 mm 等比换算（≈6.2 pt / 0.77），因此采用其**密度**而非绝对字号。

## 7. 硬件/环境备忘

- master（联网）：Foldseek 搜索、PDB 下载、数据分析、绘图
- gpu3（GPU）：RAPiDock 对接；conda 环境 `RAPiDock`；路径映射 `/dataN/zhanghaohong` → gpu3 `/data(N+3)/zhanghaohong`（home 不变）
- gpu3 有 8× A100-80GB；本次竞赛对比用 4–5 卡并行，约 20 分钟完成 327 个复合物
- SSH 用 paramiko（本机无 sshpass）
- 原始对接结果在 `peptide_docking/RAPiDock/results/`（不在本包内）
# NK-RCPP
