#!/bin/bash
# ===========================================================================
# NK1R 和 NK3R Foldseek 搜索 + 解析 + 下载 PDB 全流程（需要联网）
#
# 输出目录不在仓库内（Foldseek 结果与下载的结构体积很大），默认写到
# $HOME/nk_rcpp_work/peptide_docking，可用环境变量覆盖:
#   PDB_DIR       受体结构目录，默认仓库内的 step2_positive_docking/input/pdbs
#   DOCKING_DIR   Foldseek / 对接工作目录
#   PEPTIDES_DIR  PepMLM 导出的候选肽 CSV 所在目录（首列为序列、含表头）
# ===========================================================================
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "${SCRIPT_DIR}")"
PDB_DIR="${PDB_DIR:-${BASE_DIR}/step2_positive_docking/input/pdbs}"
DOCKING_DIR="${DOCKING_DIR:-${NK_RCPP_WORK:-$HOME/nk_rcpp_work}/peptide_docking}"
PEPTIDES_DIR="${PEPTIDES_DIR:-${DOCKING_DIR}}"

echo "=========================================="
echo "Step 1: Foldseek Search for NK1R"
echo "=========================================="
python ${SCRIPT_DIR}/foldseek_search.py \
    --pdb ${PDB_DIR}/NK1R.pdb \
    --output ${DOCKING_DIR}/nk1r_foldseek \
    --databases "pdb100,afdb50,afdb-swissprot,afdb-proteome"

echo ""
echo "=========================================="
echo "Step 2: Foldseek Search for NK3R"
echo "=========================================="
python ${SCRIPT_DIR}/foldseek_search.py \
    --pdb ${PDB_DIR}/NK3R.pdb \
    --output ${DOCKING_DIR}/nk3r_foldseek \
    --databases "pdb100,afdb50,afdb-swissprot,afdb-proteome"

echo ""
echo "=========================================="
echo "Step 3: Parse NK1R Results & Download PDBs"
echo "=========================================="
python ${SCRIPT_DIR}/parse_and_download.py \
    --foldseek_dir ${DOCKING_DIR}/nk1r_foldseek \
    --output_dir ${DOCKING_DIR}/nk1r_foldseek \
    --cdhit_identity 0.5

echo ""
echo "=========================================="
echo "Step 4: Parse NK3R Results & Download PDBs"
echo "=========================================="
python ${SCRIPT_DIR}/parse_and_download.py \
    --foldseek_dir ${DOCKING_DIR}/nk3r_foldseek \
    --output_dir ${DOCKING_DIR}/nk3r_foldseek \
    --cdhit_identity 0.5

echo ""
echo "=========================================="
echo "Step 5: Create NK1R Negative Docking Inputs"
echo "=========================================="
python ${SCRIPT_DIR}/create_neg_input.py \
    --receptor NK1R \
    --foldseek_dir ${DOCKING_DIR}/nk1r_foldseek/pdb_substructures \
    --peptides_csv ${PEPTIDES_DIR}/nk1r_scripts/nk1r.csv \
    --known_peptide "RPKPQQFFGLM" \
    --output_dir ${DOCKING_DIR}/RAPiDock

echo ""
echo "=========================================="
echo "Step 6: Create NK3R Negative Docking Inputs"
echo "=========================================="
python ${SCRIPT_DIR}/create_neg_input.py \
    --receptor NK3R \
    --foldseek_dir ${DOCKING_DIR}/nk3r_foldseek/pdb_substructures \
    --peptides_csv ${PEPTIDES_DIR}/nk3r_scripts/nk3r.csv \
    --known_peptide "DMHDFFVGLM" \
    --output_dir ${DOCKING_DIR}/RAPiDock

echo ""
echo "=========================================="
echo "ALL DONE!"
echo "=========================================="
