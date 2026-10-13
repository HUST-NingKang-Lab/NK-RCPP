#!/bin/bash
# 需在装好 RAPiDock 的机器上运行（需要 GPU）。可用环境变量覆盖:
#   RAPIDOCK_DIR  RAPiDock 安装目录（默认 $HOME/RAPiDock）
#   CONDA_ENV     conda 环境名（默认 RAPiDock）
RAPIDOCK_DIR="${RAPIDOCK_DIR:-$HOME/RAPiDock}"
CONDA_ENV="${CONDA_ENV:-RAPiDock}"
cd "${RAPIDOCK_DIR}" || exit 1
source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh
conda activate "${CONDA_ENV}" 2>/dev/null
export OMP_NUM_THREADS=4
# comp01 negative on GPU 1
( export CUDA_VISIBLE_DEVICES=1
  echo "[G1] comp01 $(date +%H:%M:%S)"
  python inference.py --config default_inference_args.yaml     --protein_peptide_csv comp_top10_neg_inputs/comp01_docking_input.csv     --output_dir results/COMP_TOP10_neg >> logs/fix_comp01.log 2>&1
  echo "[G1] comp01 DONE $(date +%H:%M:%S)" ) &
# missing known on GPU 2
( export CUDA_VISIBLE_DEVICES=2
  echo "[G2] known_missing $(date +%H:%M:%S)"
  python inference.py --config default_inference_args.yaml     --protein_peptide_csv comp_top10_neg_inputs/known_missing_docking_input.csv     --output_dir results/COMP_TOP10_neg >> logs/fix_known.log 2>&1
  echo "[G2] known_missing DONE $(date +%H:%M:%S)" ) &
wait
echo ALL_FIX_DONE
