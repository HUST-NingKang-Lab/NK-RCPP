#!/usr/bin/env python3
"""
解析 Foldseek pdb100 结果，CD-HIT 去冗余，下载 PDB 结构文件
在 MASTER 节点运行（需要联网下载 PDB）

用法: python parse_and_download.py --foldseek_dir <foldseek_results> --output_dir <output> [--identity 0.5]
"""
import pandas as pd
import os
import sys
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path


def parse_m8(file_path: str) -> pd.DataFrame:
    """Parse a Foldseek m8 file into a pandas DataFrame."""
    columns = [
        "query_id", "target_id", "percent_identity", "alignment_length",
        "mismatches", "gap_opens", "query_start", "query_end",
        "target_start", "target_end", "evalue", "bit_score",
        "probability", "query_length", "target_length", "query_alignment",
        "target_alignment", "target_confidence_array", "target_sequence",
        "target_taxid", "target_taxname"
    ]
    df = pd.read_csv(file_path, sep="\t", header=None, names=columns)
    return df


def main():
    parser = argparse.ArgumentParser(description="Parse Foldseek results and download PDBs")
    parser.add_argument("--foldseek_dir", required=True, help="Directory with Foldseek m8 files")
    parser.add_argument("--output_dir", required=True, help="Output directory for PDB files and substructures")
    parser.add_argument("--cdhit_identity", default="0.5", help="CD-HIT identity threshold (default: 0.5)")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    pdb_files_dir = output_dir / "pdb_files"
    substructures_dir = output_dir / "pdb_substructures"
    pdb_files_dir.mkdir(parents=True, exist_ok=True)
    substructures_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Parse pdb100 m8 file
    m8_path = os.path.join(args.foldseek_dir, "alis_pdb100.m8")
    if not os.path.exists(m8_path):
        sys.exit(f"[ERROR] {m8_path} not found")

    print(f"[INFO] Parsing {m8_path}...")
    pdb = parse_m8(m8_path)
    print(f"[INFO] Found {len(pdb)} total hits")

    # Step 2: Create FASTA from target sequences
    fasta_path = str(output_dir / "alis_pdb100.fasta")
    with open(fasta_path, "w") as fasta_file:
        for index, row in pdb.iterrows():
            seq_record = SeqRecord(
                seq=row["target_sequence"],
                id=row["target_id"],
                description=""
            )
            SeqIO.write(seq_record, fasta_file, "fasta")
    print(f"[INFO] Wrote FASTA to {fasta_path}")

    # Step 3: Run CD-HIT or fallback to Python clustering
    cdhit_fasta = str(output_dir / "alis_pdb100_cdhit50.fasta")
    cdhit_clstr = str(output_dir / "alis_pdb100_cdhit50.fasta.clstr")

    # Try CD-HIT first, fall back to Python clusterer
    import shutil
    if shutil.which("cd-hit"):
        cmd = f"cd-hit -i {fasta_path} -o {cdhit_fasta} -c {args.cdhit_identity} -n 3 -d 0"
        print(f"[INFO] Running CD-HIT: {cmd}")
        ret = os.system(cmd)
        if ret != 0:
            print(f"[WARN] CD-HIT failed, falling back to Python clusterer")
            use_python_clusterer = True
        else:
            use_python_clusterer = False
    else:
        print("[INFO] cd-hit not found, using Python clusterer")
        use_python_clusterer = True

    if use_python_clusterer:
        # Use Python k-mer based clustering
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
        py_cmd = f"python3 {os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cluster_sequences.py')} --input {fasta_path} --output {cdhit_fasta} --threshold {args.cdhit_identity}"
        print(f"[INFO] Running Python clusterer: {py_cmd}")
        ret = os.system(py_cmd)
        if ret != 0:
            sys.exit(f"[ERROR] Python clusterer failed with code {ret}")

    # Step 4: Parse CD-HIT output and extract unique PDB IDs
    cdhit_records = list(SeqIO.parse(cdhit_fasta, "fasta"))
    pdb_ids_raw = [record.id.split('-')[0] for record in cdhit_records]
    pdb_ids = list(set(pdb_ids_raw))
    print(f"[INFO] Unique PDB IDs after CD-HIT: {len(pdb_ids)}")

    # Step 5: Parse chain info from Foldseek target IDs
    # Format: pdb_id-assembly_chain, e.g., "3oe6-assembly1.cif.gz_A"
    chain_map = {}
    for record in cdhit_records:
        target_id = record.id
        parts = target_id.split('-')
        pdb_base = parts[0]
        if '_' in target_id:
            chain = target_id.split('_')[-1][0]  # First char after last underscore
            chain_map[pdb_base] = chain

    print(f"[INFO] Chain map: {len(chain_map)} entries")

    # Step 6: Download PDB files from RCSB
    print(f"[INFO] Downloading PDB files to {pdb_files_dir}...")
    downloaded = 0
    for pdb_id in pdb_ids:
        pdb_file_path = pdb_files_dir / f"{pdb_id}.pdb"
        if not pdb_file_path.exists():
            cmd = f"wget -q https://files.rcsb.org/download/{pdb_id}.pdb -O {pdb_file_path}"
            ret = os.system(cmd)
            if ret == 0:
                downloaded += 1
                print(f"  Downloaded {pdb_id}.pdb")
            else:
                print(f"  [WARN] Failed to download {pdb_id}.pdb")
        else:
            print(f"  [SKIP] {pdb_id}.pdb already exists")
    print(f"[INFO] Downloaded {downloaded} new PDB files")

    # Step 7: Extract substructures (chains)
    print(f"[INFO] Extracting chain substructures...")
    from Bio.PDB import PDBParser, PDBIO
    pdb_parser = PDBParser(QUIET=True)
    pdb_io = PDBIO()

    # Also extract chains for PDBs not in chain_map (use first chain)
    extracted = 0
    for pdb_file in sorted(pdb_files_dir.glob("*.pdb")):
        complex_id = pdb_file.stem
        try:
            structure = pdb_parser.get_structure(complex_id, str(pdb_file))
        except Exception as e:
            print(f"  [WARN] Failed to parse {pdb_file.name}: {e}")
            continue

        for model in structure:
            chain_ids = [chain.id for chain in model]

            if complex_id in chain_map and chain_map[complex_id] in chain_ids:
                target_chain = chain_map[complex_id]
            else:
                target_chain = chain_ids[0]  # Use first chain as fallback

            output_path = substructures_dir / f"{complex_id}_{target_chain}.pdb"
            if not output_path.exists():
                # Write only the target chain
                class ChainSelector:
                    def __init__(self, chain_id):
                        self.chain_id = chain_id
                    def accept_chain(self, chain):
                        return chain.id == self.chain_id
                    def accept_model(self, model):
                        return True
                    def accept_residue(self, residue):
                        return True
                    def accept_atom(self, atom):
                        return True

                pdb_io.set_structure(structure)
                pdb_io.save(str(output_path), ChainSelector(target_chain))
                extracted += 1
                print(f"  Extracted {output_path.name}")

    print(f"[INFO] Extracted {extracted} substructures")
    print(f"[INFO] Done! Results in {output_dir}")


if __name__ == "__main__":
    main()
