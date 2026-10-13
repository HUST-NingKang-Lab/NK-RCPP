#!/bin/bash
# ===========================================================================
# 竞赛 Top10 肽的 RAPiDock 对接 (NK2R)
# 1. 正向: Top10 vs NK2R      -> results/COMP_TOP10_NK2R/{comp01..comp10}/
# 2. 负向: Top10 + 已知肽 vs 27 个脱靶蛋白 -> results/COMP_TOP10_neg/
# 可重复执行 (自动跳过已完成批次)。需要 GPU。
#
# 环境变量:
#   RAPIDOCK_DIR  RAPiDock 安装目录，默认 $HOME/RAPiDock
#   CONDA_ENV     conda 环境名，默认 RAPiDock
# ===========================================================================
set -u

# 需在装好 RAPiDock 的机器上运行（需要 GPU）。可用环境变量覆盖:
#   RAPIDOCK_DIR  RAPiDock 安装目录（默认 $HOME/RAPiDock）
#   CONDA_ENV     conda 环境名（默认 RAPiDock）
RAPIDOCK_DIR="${RAPIDOCK_DIR:-$HOME/RAPiDock}"
CONDA_ENV="${CONDA_ENV:-RAPiDock}"
cd "${RAPIDOCK_DIR}" || exit 1
source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh
conda activate "${CONDA_ENV}" 2>/dev/null

LOG_DIR="${RAPIDOCK_DIR}/logs"
mkdir -p "${LOG_DIR}"

echo "=========== START $(date) ==========="

# ---------- 1. Positive docking ----------
echo "--- Positive docking: Top10 vs NK2R ---"
mkdir -p results/COMP_TOP10_NK2R
for f in comp_top10_pos_inputs/*.csv; do
    name=$(basename "$f" _docking_input.csv)
    if [ -f "results/COMP_TOP10_NK2R/${name}/ref2015_score.csv" ]; then
        echo "[SKIP] ${name} (done)"; continue
    fi
    echo "[RUN] ${name} $(date)"
    python inference.py --config default_inference_args.yaml \
        --protein_peptide_csv "$f" \
        --output_dir results/COMP_TOP10_NK2R >> "${LOG_DIR}/comp_pos_${name}.log" 2>&1
    echo "[DONE] ${name} $(date)"
done

# ---------- 2. Negative docking ----------
echo "--- Negative docking: Top10 + known vs 27 off-targets ---"
mkdir -p results/COMP_TOP10_neg
for f in comp_top10_neg_inputs/*.csv; do
    name=$(basename "$f" _docking_input.csv)
    n_done=$(find "results/COMP_TOP10_neg/${name}" -name ref2015_score.csv 2>/dev/null | wc -l)
    if [ "${n_done}" -ge 27 ]; then
        echo "[SKIP] ${name} (${n_done}/27 done)"; continue
    fi
    echo "[RUN] ${name} (${n_done}/27 already done) $(date)"
    python inference.py --config default_inference_args.yaml \
        --protein_peptide_csv "$f" \
        --output_dir results/COMP_TOP10_neg >> "${LOG_DIR}/comp_neg_${name}.log" 2>&1
    echo "[DONE] ${name} $(date)"
done

echo "=========== ALL DONE $(date) ==========="
