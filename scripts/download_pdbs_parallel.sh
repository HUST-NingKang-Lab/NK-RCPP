#!/bin/bash
# 并行下载 PDB 文件（8线程）
# 用法: bash download_pdbs_parallel.sh <pdb_ids_file> <output_dir>
set -e

PDB_LIST=$1
OUTPUT_DIR=$2
THREADS=8

download_one() {
    pdb_id=$1
    output_dir=$2
    pdb_file="${output_dir}/${pdb_id}.pdb"
    if [ ! -f "$pdb_file" ]; then
        wget -q --timeout=30 -O "$pdb_file" "https://files.rcsb.org/download/${pdb_id}.pdb" 2>/dev/null && echo "  OK: $pdb_id" || echo "  FAIL: $pdb_id"
    fi
}
export -f download_one

echo "[INFO] Downloading $(wc -l < $PDB_LIST) PDBs to $OUTPUT_DIR (${THREADS} threads)..."
cat "$PDB_LIST" | xargs -P $THREADS -I {} bash -c 'download_one "$@"' _ {} "$OUTPUT_DIR"
echo "[INFO] Done!"
