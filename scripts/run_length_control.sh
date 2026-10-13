#!/bin/bash
# 需在装好 RAPiDock 的机器上运行（需要 GPU）。可用环境变量覆盖:
#   RAPIDOCK_DIR  RAPiDock 安装目录（默认 $HOME/RAPiDock）
#   CONDA_ENV     conda 环境名（默认 RAPiDock）
RAPIDOCK_DIR="${RAPIDOCK_DIR:-$HOME/RAPiDock}"
CONDA_ENV="${CONDA_ENV:-RAPiDock}"
cd "${RAPIDOCK_DIR}" || exit 1
source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/miniconda3/etc/profile.d/conda.sh
conda activate "${CONDA_ENV}" 2>/dev/null
export CUDA_VISIBLE_DEVICES=4; export OMP_NUM_THREADS=4
mkdir -p results/LENGTH_CONTROL
for f in length_control_inputs/*.csv; do
  n=$(basename $f _docking_input.csv)
  [ -f results/LENGTH_CONTROL/$n/ref2015_score.csv ] && { echo "[SKIP] $n"; continue; }
  echo "[RUN] $n $(date +%H:%M:%S)"
  python inference.py --config default_inference_args.yaml --protein_peptide_csv $f --output_dir results/LENGTH_CONTROL >> logs/lenctl_$n.log 2>&1
  echo "[DONE] $n $(date +%H:%M:%S)"
done
echo "LENGTH CONTROL FINISHED"
