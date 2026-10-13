#!/usr/bin/env python3
"""
为 NK1R/NK3R 创建负向对接输入 CSV 文件
用法: python create_neg_input.py --receptor {NK1R,NK3R} --foldseek_dir <path> --peptides_csv <path> --output <path>
"""
import argparse
import pandas as pd
import os
from pathlib import Path


def list_pdb_files(pdb_dir: Path) -> list:
    """List all .pdb files in a directory."""
    if not pdb_dir.exists():
        print(f"[WARN] PDB directory not found: {pdb_dir}")
        return []
    files = sorted([p for p in pdb_dir.glob("*.pdb") if p.is_file()])
    if not files:
        print(f"[WARN] No .pdb files found in {pdb_dir}")
    return files


def check_pdb_file(pdb_path: Path) -> bool:
    """Check if PDB file is valid (no TER in middle)."""
    with open(pdb_path, 'r') as f:
        lines = f.readlines()
    ter_indices = [i for i, line in enumerate(lines) if line.startswith('TER')]
    if len(ter_indices) > 1:
        print(f"  [SKIP] Multiple TER lines in {pdb_path.name}")
        return False
    if len(ter_indices) == 1 and ter_indices[0] < len(lines) - 3:
        print(f"  [SKIP] TER line in middle of {pdb_path.name}")
        return False
    return True


def main():
    parser = argparse.ArgumentParser(description="Create negative docking input CSVs")
    parser.add_argument("--receptor", required=True, choices=["NK1R", "NK3R"],
                        help="Receptor name")
    parser.add_argument("--foldseek_dir", required=True,
                        help="Path to foldseek pdb_substructures directory")
    parser.add_argument("--peptides_csv", required=True,
                        help="Path to candidate peptides CSV (with Pseudo Perplexity)")
    parser.add_argument("--known_peptide", required=True,
                        help="Known peptide sequence (control)")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory for docking input CSVs")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load peptides
    peptides_df = pd.read_csv(args.peptides_csv, header=None)
    peptides = peptides_df[0].tolist()[1:]  # Skip header
    print(f"[INFO] Loaded {len(peptides)} candidate peptides for {args.receptor}")

    # Get valid PDB substructures
    pdb_dir = Path(args.foldseek_dir)
    pdb_files = list_pdb_files(pdb_dir)
    print(f"[INFO] Found {len(pdb_files)} PDB substructures")

    pdb_files = [p for p in pdb_files if check_pdb_file(p)]
    print(f"[INFO] {len(pdb_files)} valid PDB files after filtering")

    if not pdb_files:
        print("[ERROR] No valid PDB files found!")
        return

    # Build dataframe: each peptide × each protein
    rows = []
    for i, pep in enumerate(peptides, start=1):
        label = f"peptide{i:02d}"
        for pdb_path in pdb_files:
            pdb_id = pdb_path.stem
            complex_name = f"{pdb_id}_{label}"
            rows.append({
                "complex_name": complex_name,
                "protein_description": str(pdb_path.absolute()),
                "peptide_description": pep,
            })

    df = pd.DataFrame(rows, columns=["complex_name", "protein_description", "peptide_description"])
    print(f"[INFO] Total docking pairs: {len(df)} ({len(peptides)} peptides × {len(pdb_files)} proteins)")

    # Also add known peptide
    known_rows = []
    for pdb_path in pdb_files:
        pdb_id = pdb_path.stem
        complex_name = f"{pdb_id}_known"
        known_rows.append({
            "complex_name": complex_name,
            "protein_description": str(pdb_path.absolute()),
            "peptide_description": args.known_peptide,
        })
    known_df = pd.DataFrame(known_rows)

    # Save combined file
    full_df = pd.concat([df, known_df], ignore_index=True)
    combined_path = output_dir / f"{args.receptor}_negative_docking_input.csv"
    full_df.to_csv(combined_path, index=False)
    print(f"[INFO] Saved combined input to {combined_path}")

    # Split by peptide into individual CSVs
    neg_inputs_dir = output_dir / f"{args.receptor.lower()}_neg_inputs"
    neg_inputs_dir.mkdir(parents=True, exist_ok=True)

    for i, pep in enumerate(peptides, start=1):
        label = f"peptide{i:02d}"
        sub_df = df[df['complex_name'].str.contains(f"_{label}$")]
        sub_path = neg_inputs_dir / f"peptide{i:02d}_docking_input.csv"
        sub_df.to_csv(sub_path, index=False)

    # Known peptide
    known_path = neg_inputs_dir / "known_peptide_docking_input.csv"
    known_df.to_csv(known_path, index=False)

    print(f"[INFO] Split into {len(peptides) + 1} individual CSVs in {neg_inputs_dir}")
    print(f"[INFO] Done!")


if __name__ == "__main__":
    main()
