#!/usr/bin/env python3
"""
从竞赛排行榜 Excel 提取 Top10 肽，生成 RAPiDock 对接输入:
1. 正向对接: Top10 vs NK2R.pdb
2. 负向对接: Top10 vs NK2R 的 27 个脱靶蛋白 (nk2r_foldseek/pdb_substructures)
输出到 peptide_docking/RAPiDock/comp_top10_* 目录
"""
import pandas as pd
import os
from pathlib import Path

# --- paths ---------------------------------------------------------------
# Repository-relative: the scripts run from anywhere inside a checkout.
# Raw RAPiDock run directories and the competition workbook are not shipped;
# set NK_RCPP_RAW / NK_RCPP_COMP_XLSX (or NK_RCPP_BASE) to point at them.
# See scripts/README.md.
BASE = Path(os.environ.get("NK_RCPP_BASE", Path(__file__).resolve().parent.parent))
ANALYSIS = BASE / "step5_analysis" / "output"
RAW = Path(os.environ.get("NK_RCPP_RAW", BASE / "raw"))
DOCKING = RAW / "peptide_docking"
RAPIDOCK = DOCKING / "RAPiDock"
XLSX = Path(os.environ.get("NK_RCPP_COMP_XLSX",
                  BASE / "step6_competition_benchmark" / "input" /
                  "2025_pep_design_comp_all_entries_with_synthesize_highlight.xlsx"))
# ---------------------------------------------------------------------------

RAPIDOCK = DOCKING / "RAPiDock"
NK2R_PDB = DOCKING / "pdbs/NK2R.pdb"
OFFTARGET_DIR = DOCKING / "nk2r_foldseek/pdb_substructures"
OUT_POS_DIR = RAPIDOCK / "comp_top10_pos_inputs"
OUT_NEG_DIR = RAPIDOCK / "comp_top10_neg_inputs"

RANK_COL = "Peptide Name"
SEQ_COL = "Amino acid Sequence"
EC50_COL = "EC50 Value (nM)"
TOP_N = 10


def check_pdb_file(pdb_path: Path) -> bool:
    """Reject PDBs with TER in the middle (same filter as our pipeline)."""
    with open(pdb_path) as f:
        lines = f.readlines()
    ter_indices = [i for i, line in enumerate(lines) if line.startswith("TER")]
    if len(ter_indices) > 1:
        return False
    if len(ter_indices) == 1 and ter_indices[0] < len(lines) - 3:
        return False
    return True


def main():
    df = pd.read_excel(XLSX, sheet_name="2025_pep_design_comp_all_entrie")
    # Leaderboard order = file row order (verified: sorted by EC50 ascending)
    top = df.head(TOP_N).copy()

    top = top.reset_index(drop=True)
    print(f"=== Competition Top {TOP_N} (leaderboard order) ===")
    for i in range(len(top)):
        code = f"COMP-{i+1:02d}"
        name = str(top.loc[i, RANK_COL])
        seq = str(top.loc[i, SEQ_COL]).strip()
        ec50 = top.loc[i, EC50_COL]
        print(f"  {code}: {name:20s} {seq:35s} EC50={ec50} nM")

    # ---- 1. Positive docking: Top10 vs NK2R ----
    OUT_POS_DIR.mkdir(parents=True, exist_ok=True)
    pos_rows = []
    for i in range(len(top)):
        code = f"comp{i+1:02d}"
        pos_rows.append({
            "complex_name": code,
            "protein_description": str(NK2R_PDB),
            "peptide_description": str(top.loc[i, SEQ_COL]).strip(),
        })
    pos_df = pd.DataFrame(pos_rows, columns=["complex_name", "protein_description", "peptide_description"])
    pos_csv = RAPIDOCK / "COMP_TOP10_NK2R_peptide.csv"
    pos_df.to_csv(pos_csv, index=False)
    print(f"\n[OK] Positive input: {pos_csv} ({len(pos_df)} peptides)")

    # Split per peptide for robust resume
    for i in range(len(top)):
        label = f"comp{i+1:02d}"
        sub = pos_df[pos_df["complex_name"] == label]
        sub.to_csv(OUT_POS_DIR / f"{label}_docking_input.csv", index=False)

    # ---- 2. Negative docking: Top10 vs off-targets ----
    OUT_NEG_DIR.mkdir(parents=True, exist_ok=True)
    off_targets = sorted(p for p in OFFTARGET_DIR.glob("*.pdb") if check_pdb_file(p))
    print(f"[INFO] Off-target proteins: {len(off_targets)} valid")
    if not off_targets:
        print("[ERROR] no valid off-target PDBs")
        return

    neg_rows = []
    for i in range(len(top)):
        label = f"comp{i+1:02d}"
        seq = str(top.loc[i, SEQ_COL]).strip()
        for pdb_path in off_targets:
            neg_rows.append({
                "complex_name": f"{pdb_path.stem}_comp{i:02d}",
                "protein_description": str(pdb_path.absolute()),
                "peptide_description": seq,
            })
            # known ligand row for paired comparison (same as our pipeline)
            neg_rows.append({
                "complex_name": f"{pdb_path.stem}_known",
                "protein_description": str(pdb_path.absolute()),
                "peptide_description": "HKTDSFVGLM",
            })
    neg_df = pd.DataFrame(neg_rows, columns=["complex_name", "protein_description", "peptide_description"]).drop_duplicates()
    neg_csv = RAPIDOCK / "COMP_TOP10_negative_docking_input.csv"
    neg_df.to_csv(neg_csv, index=False)
    print(f"[OK] Negative input: {neg_csv} ({len(neg_df)} pairs)")

    # split per peptide
    for i in range(1, TOP_N + 1):
        label = f"comp{i:02d}"
        sub = neg_df[neg_df["complex_name"].str.contains(f"_{label}$")]
        sub.to_csv(OUT_NEG_DIR / f"{label}_docking_input.csv", index=False)
    # known ligand batch
    known_sub = neg_df[neg_df["complex_name"].str.endswith("_known")]
    known_sub.to_csv(OUT_NEG_DIR / "known_peptide_docking_input.csv", index=False)
    print(f"[OK] Split into {TOP_N} peptide batches + 1 known batch in {OUT_NEG_DIR}")

    # save metadata for later comparison
    meta = pd.DataFrame({
        "code": [f"COMP-{i:02d}" for i in range(1, TOP_N + 1)],
        "peptide_name": top[RANK_COL].astype(str).tolist(),
        "sequence": top[SEQ_COL].astype(str).str.strip().tolist(),
        "EC50_nM": top[EC50_COL].tolist(),
        "rank": list(range(1, TOP_N + 1)),
    })
    meta.to_csv(ANALYSIS / "competition_top10.csv", index=False)
    print(f"[OK] Metadata saved to analysis/competition_top10.csv")


if __name__ == "__main__":
    main()
