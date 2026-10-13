#!/usr/bin/env python3
"""
提取 PDB 子链结构
用法: python extract_chains.py --foldseek_dir <path> --output_dir <path>
"""
import argparse
import os
from pathlib import Path
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--foldseek_dir", required=True, help="Foldseek results directory")
    parser.add_argument("--output_dir", required=True, help="Output directory for substructures")
    args = parser.parse_args()

    foldseek_dir = Path(args.foldseek_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Read CD-HIT results to get chain info
    cdhit_fasta = foldseek_dir / "alis_pdb100_cdhit50.fasta"
    records = list(SeqIO.parse(str(cdhit_fasta), "fasta"))

    # Build chain map: pdb_id -> chain letter
    chain_map = {}
    for record in records:
        target_id = record.id
        pdb_base = target_id.split('-')[0]
        if '_' in target_id:
            chain = target_id.split('_')[-1][0]
            chain_map[pdb_base] = chain
    print(f"[INFO] Chain map: {len(chain_map)} entries")

    # Parse and extract
    pdb_files_dir = foldseek_dir / "pdb_files"
    pdb_parser = PDBParser(QUIET=True)
    pdb_io = PDBIO()

    extracted = 0
    failed = 0
    skipped = 0

    for pdb_file in sorted(pdb_files_dir.glob("*.pdb")):
        complex_id = pdb_file.stem
        output_path = output_dir / f"{complex_id}_{chain_map.get(complex_id, 'A')}.pdb"

        if output_path.exists():
            skipped += 1
            continue

        try:
            structure = pdb_parser.get_structure(complex_id, str(pdb_file))
        except Exception as e:
            print(f"  [FAIL] Parse {pdb_file.name}: {e}")
            failed += 1
            continue

        chain_ids = [chain.id for chain in structure[0]]
        target_chain = chain_map.get(complex_id, chain_ids[0])

        if target_chain not in chain_ids:
            target_chain = chain_ids[0]
            print(f"  [WARN] Chain {chain_map.get(complex_id, '?')} not in {complex_id}, using {target_chain}")

        try:
            pdb_io.set_structure(structure)
            pdb_io.save(str(output_path))
            extracted += 1
        except Exception as e:
            print(f"  [FAIL] Extract {pdb_file.name}: {e}")
            failed += 1

    print(f"[INFO] Extracted: {extracted}, Skipped: {skipped}, Failed: {failed}")


if __name__ == "__main__":
    main()
