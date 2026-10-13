#!/bin/bash
# 竞赛 Top10 对接 worker — 指定 GPU + 批次列表
# 用法: bash comp_top10_worker.sh <GPU_ID> <phase> <batch1> <batch2> ...
#   phase: pos | neg
set -u

GPU_ID=$1; shift
PHASE=$1; shift
BATCHES="$@"

# 需在装好 RAPiDock 的机器上运行（需要 GPU）。可用环境变量覆盖:
#   RAPIDOCK_DIR  RAPiDock 安装目录（默认 $HOME/RAPiDock）
#   CONDA_ENV     conda 环境名（默认 RAPiDock）
RAPIDOCK_DIR="${RAPIDOCK_DIR:-$HOME/RAPiDock}"
CONDA_ENV="${CONDA_ENV:-RAPiDock}"
cd "${RAPIDOCK_DIR}" || exit 1
source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh
conda activate "${CONDA_ENV}" 2>/dev/null
export CUDA_VISIBLE_DEVICES=${GPU_ID}
export OMP_NUM_THREADS=4

mkdir -p logs

if [ "${PHASE}" = "pos" ]; then
    IN_DIR="comp_top10_pos_inputs"; OUT_DIR="results/COMP_TOP10_NK2R"; N_EXPECT=1
else
    IN_DIR="comp_top10_neg_inputs"; OUT_DIR="results/COMP_TOP10_neg"; N_EXPECT=27
fi
mkdir -p "${OUT_DIR}"

for name in ${BATCHES}; do
    f="${IN_DIR}/${name}_docking_input.csv"
    [ -f "$f" ] || { echo "[GPU${GPU_ID}] [MISS] $f"; continue; }
    n_done=$(find "${OUT_DIR}/${name}" -name ref2015_score.csv 2>/dev/null | wc -l)
    if [ "${n_done}" -ge "${N_EXPECT}" ]; then
        echo "[GPU${GPU_ID}] [SKIP] ${name} (${n_done}/${N_EXPECT})"; continue
    fi
    echo "[GPU${GPU_ID}] [RUN] ${name} $(date +%H:%M:%S)"
    python inference.py --config default_inference_args.yaml \
        --protein_peptide_csv "$f" \
        --output_dir "${OUT_DIR}" >> "logs/gpu${GPU_ID}_${name}.log" 2>&1
    echo "[GPU${GPU_ID}] [DONE] ${name} $(date +%H:%M:%S)"
done
echo "[GPU${GPU_ID}] WORKER FINISHED"
