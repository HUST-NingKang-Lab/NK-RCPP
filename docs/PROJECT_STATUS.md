# 项目状态与数据来源说明

> 本文件记录**数据来源与运行环境**，便于复现；结论与数字一律以稿件
> *NK-RCPP: a receptor-conditioned computational workflow for prioritizing neurokinin-receptor peptide
> candidates* 为准（见 [`MANUSCRIPT_GUIDE.md`](MANUSCRIPT_GUIDE.md) 速查表）。
> ⚠️ 早期阶段的临时数字不要在本仓库中继续引用。

## 流程

```
PepMLM 受体条件化生成 → RAPiDock 靶内对接（vs 受体）
                       → Foldseek 结构近邻检索
                       → RAPiDock 结构近邻对接（候选与对照同一近邻面板）
                       → 确定性优先级规则（primary / strict）
```

## 受体匹配对照

| 受体 | 对照 | 序列 |
|------|------|------|
| NK1R | Substance P | `RPKPQQFFGLM` |
| NK2R | NKA 衍生肽 | `HKTDSFVGLM` |
| NK3R | Neurokinin B | `DMHDFFVGLM` |

## 运行环境与路径

* Conda 环境：`RAPiDock`；脚本通过 `RAPIDOCK_DIR` 环境变量定位 RAPiDock 安装目录（默认 `$HOME/RAPiDock`）
* 对接机器**无外网**：Foldseek 走 Web API（在有网的机器上完成），PDB 下载由
  `scripts/download_pdbs_parallel.sh` 预先完成
* 中间产物（Foldseek 结果、下载的结构、RAPiDock 运行目录）体积大，不入库；用 `NK_RCPP_RAW` 等环境变量指向本地副本
* 凭证不写入本仓库；访问远程机器请使用密钥或个人配置

## 对接设置

* 每条肽 10 个构象（`N = 10`），16 步扩散（`inference_steps` / `actual_steps = 16`）
* 打分函数 `ref2015`，本地权重 `rapidock_local.pt`，`conformation_partial = 1:1:1`，`batch_size = 4`
* 归档聚合取可得的**最低** `ref2015score`；下游分析统一使用 `rank1`

## 完成状态

| 步骤 | NK1R | NK2R | NK3R | 对照实验 |
|------|------|------|------|----------|
| PepMLM 生成 | ✅ 32 | ✅ 32 | ✅ 32 | — |
| 靶内对接 | ✅ | ✅ | ✅ | ✅ COMP Top10（10）+ 长度对照（20） |
| Foldseek 检索 | ✅ | ✅ | ✅ | — |
| 结构近邻对接 | ✅ 66 × 33 | ✅ 27 × 19 | ✅ 60 × 33 | ✅ 27 × 11 |
| 审计聚合 | ✅ | ⚠️ 18/32 可评估 | ✅ | ✅ 探索性 |

## ⚠️ NK2R 数据缺口（稿件中的缺失性边界）

14 条 NK2R 候选没有结构近邻聚合值（原始 `valid_pep` 子集只选了 18 条候选跑结构近邻对接）：

```
pep2, pep5, pep6, pep8, pep12, pep15, pep16, pep17, pep18, pep20, pep21, pep25, pep29, pep31
```

其中 10 条 Delta-on < 0。由此得到稿件中的无约束边界：**primary 10–20、strict 2–12**。
该边界只量化计数对不可得证据的依赖，**不预测**缺失结果。

补跑方式：由 `scripts/create_neg_input.py` 生成缺失批次的 `nk2r_neg_inputs/*.csv`，在 GPU 节点运行结构近邻
对接，再用 `scripts/harvest_docking_scores.py` 与 `scripts/pairwise_selectivity.py` 重新汇总。

## 本仓库相对被审计软件包多出的内容

被审计的软件包只有候选级聚合值。本仓库保留了：

* **逐复合物结构近邻打分**（[`../step4_negative_docking/output/negative_scores_long.csv`](../step4_negative_docking/output/)）
  → 可重算全部 4,518 次配对比较 / 82 条可评估候选（`scripts/paired_neighbour_deltas.py`，与归档聚合零分歧）
* **Foldseek 命中表**（[`../step3_foldseek_search/output/`](../step3_foldseek_search/)）→ 面板成员身份可查
* **竞赛与长度对照的原始打分**（[`../step6_competition_benchmark/output/`](../step6_competition_benchmark/)）

仍不可得：RAPiDock 原始构象目录、原始生成输出（含 NK2R 伪困惑度）、COMP/长度对照战役的逐对行。

## 待办

1. 若需补全 NK2R 的联合规则留存统计：补跑上述 14 条候选的结构近邻对接
2. 前瞻性验证：等规模、受体特异的 testing budget；strict 的 2 条 NK2R 候选作为预设亚组
