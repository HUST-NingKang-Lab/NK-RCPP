#!/bin/bash
# ===========================================================================
# NK1R 和 NK3R 负向对接脚本
# 需要 GPU。
#
# 用法:
#   1. 先完成 Foldseek 搜索 + PDB 下载 + 对接输入 CSV 生成
#      （见 scripts/run_full_pipeline_master.sh）
#   2. 把 nk1r_neg_inputs/ 与 nk3r_neg_inputs/ 放到 $RAPIDOCK_DIR 下
#   3. bash scripts/run_docking_gpu.sh
#
# 环境变量:
#   RAPIDOCK_DIR  RAPiDock 安装目录（默认 $HOME/RAPiDock）
#   CONDA_ENV     conda 环境名（默认 RAPiDock）
# ===========================================================================
set -e

RAPIDOCK_DIR="${RAPIDOCK_DIR:-$HOME/RAPiDock}"
CONDA_ENV="${CONDA_ENV:-RAPiDock}"

# 激活 conda 环境
source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh
conda activate "${CONDA_ENV}"

cd "${RAPIDOCK_DIR}" || exit 1

echo "=============================================="
echo "NK1R Negative Docking"
echo "=============================================="
mkdir -p results/NK1R_neg

# 处理每个 peptide
for f in nk1r_neg_inputs/*.csv; do
    name=$(basename "$f" .csv)
    output_dir="results/NK1R_neg/${name}"

    if [ -d "$output_dir" ]; then
        # 检查是否有 ref2015_score.csv，有则跳过
        completed=$(find "$output_dir" -name "ref2015_score.csv" 2>/dev/null | wc -l)
        if [ "$completed" -gt 0 ]; then
            echo "[SKIP] $name already completed ($completed sub-complexes)"
            continue
        fi
    fi

    echo "[RUNNING] $name..."
    python inference.py --config default_inference_args.yaml \
        --protein_peptide_csv "$f" \
        --output_dir "$output_dir"
    echo "[DONE] $name"
done

echo ""
echo "=============================================="
echo "NK3R Negative Docking"
echo "=============================================="
mkdir -p results/NK3R_neg

for f in nk3r_neg_inputs/*.csv; do
    name=$(basename "$f" .csv)
    output_dir="results/NK3R_neg/${name}"

    if [ -d "$output_dir" ]; then
        completed=$(find "$output_dir" -name "ref2015_score.csv" 2>/dev/null | wc -l)
        if [ "$completed" -gt 0 ]; then
            echo "[SKIP] $name already completed ($completed sub-complexes)"
            continue
        fi
    fi

    echo "[RUNNING] $name..."
    python inference.py --config default_inference_args.yaml \
        --protein_peptide_csv "$f" \
        --output_dir "$output_dir"
    echo "[DONE] $name"
done

echo ""
echo "=============================================="
echo "ALL DONE!"
echo "=============================================="
