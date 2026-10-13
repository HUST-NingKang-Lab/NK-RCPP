#!/usr/bin/env python3
"""
生成论文附表:
- Table S1: 代号 → 序列映射 + PepMLM 信息
- Table S2: 全部候选肽完整对接数据 (正/负/配对选择性)
- Table S3: 脱靶蛋白身份表 (PDB, 名称, 蛋白类型)
"""
import pandas as pd
import numpy as np
import glob
import os
from pathlib import Path

# --- paths ---------------------------------------------------------------
# Repository-relative: the scripts run from anywhere inside a checkout.
# Raw RAPiDock run directories and the competition workbook are not shipped;
# set NK_RCPP_RAW / NK_RCPP_COMP_XLSX (or NK_RCPP_BASE) to point at them.
# See scripts/README.md.
BASE = Path(os.environ.get("NK_RCPP_BASE", Path(__file__).resolve().parent.parent))
ANALYSIS = BASE / "step5_analysis" / "output"
FIGDIR = BASE / "figures"
TABLES = FIGDIR / "tables"
RAW = Path(os.environ.get("NK_RCPP_RAW", BASE / "raw"))
DOCKING = RAW / "peptide_docking"
RAPIDOCK = DOCKING / "RAPiDock"
# ---------------------------------------------------------------------------

OUT = TABLES
OUT.mkdir(parents=True, exist_ok=True)

KNOWN = {"NK1R": "RPKPQQFFGLM", "NK2R": "HKTDSFVGLM", "NK3R": "DMHDFFVGLM"}

# ---------------- Table S1: code mapping ----------------
def tableS1():
    rows = []
    for rec in ["NK1R", "NK2R", "NK3R"]:
        df = pd.read_csv(ANALYSIS / f"{rec}_summary.csv")
        for _, r in df.iterrows():
            pid = r["peptide_id"]
            seq = str(r["sequence"]).strip()
            if pid in ["Substance_P", "pep_known", "Neurokinin_B"]:
                code, ligand = f"{rec}-K", "Yes"
            else:
                code, ligand = f"{rec}-{int(pid.replace('pep','')):02d}", "No"
            rows.append({"Code": code, "Receptor": rec, "Ligand": ligand,
                         "PepMLM_ID": pid, "Sequence": seq})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "TableS1_peptide_codes.csv", index=False)
    print(f"[OK] TableS1: {len(df)} rows")

# ---------------- Table S2: full docking data ----------------
def tableS2():
    rows = []
    for rec in ["NK1R", "NK2R", "NK3R"]:
        rank = pd.read_csv(ANALYSIS / f"{rec}_final_ranking.csv")
        for _, r in rank.iterrows():
            pid = r["peptide_id"]
            seq = str(r["sequence"]).strip()
            if pid in ["Substance_P", "pep_known", "Neurokinin_B"]:
                code = f"{rec}-K"
            else:
                code = f"{rec}-{int(pid.replace('pep','')):02d}"
            rows.append({
                "Code": code, "Receptor": rec, "Sequence": seq,
                "Pos_ref2015": round(r["pos_score"], 2) if pd.notna(r["pos_score"]) else np.nan,
                "Neg_mean": round(r["neg_mean"], 2) if pd.notna(r["neg_mean"]) else np.nan,
                "Neg_median": round(r["neg_median"], 2) if pd.notna(r["neg_median"]) else np.nan,
                "N_offtarget": int(r["n_neg"]) if pd.notna(r["n_neg"]) else 0,
                "DeltaPos": round(r["ΔPos"], 2) if pd.notna(r.get("ΔPos", np.nan)) else np.nan,
                "DeltaNeg_mean": round(r["ΔNeg"], 2) if pd.notna(r.get("ΔNeg", np.nan)) else np.nan,
                "Frac_safer": round(r["frac_safer"], 3) if pd.notna(r.get("frac_safer", np.nan)) else np.nan,
                "Better_binding": r.get("better_binding", np.nan),
                "Safer": r.get("safer", np.nan),
                "Both": r.get("both", np.nan),
            })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "TableS2_docking_data.csv", index=False)
    print(f"[OK] TableS2: {len(df)} rows")

# ---------------- Table S3: off-target proteins ----------------
def tableS3():
    rows = []
    for rec in ["NK1R", "NK2R", "NK3R"]:
        # Proteins actually used in negative docking = unique from neg input CSVs
        neg_inputs = RAPIDOCK / f"{rec.lower()}_neg_inputs"
        proteins = set()
        # use peptide01 for the full protein list (all peptides share the same proteins)
        f1 = neg_inputs / "peptide01_docking_input.csv"
        if f1.exists():
            df1 = pd.read_csv(f1)
            proteins.update(c.split("_peptide")[0] for c in df1["complex_name"])
        fk = neg_inputs / "known_peptide_docking_input.csv"
        if fk.exists():
            dfk = pd.read_csv(fk)
            proteins.update(c.split("_known")[0] for c in dfk["complex_name"])

        sub_dir = DOCKING / f"{rec.lower()}_foldseek/pdb_substructures"
        orig_dir = DOCKING / f"{rec.lower()}_foldseek/pdb_files"
        for pdb_chain in sorted(proteins):
            pdb_id = pdb_chain.split("_")[0]
            # TITLE lives in the original full PDB, not the extracted chain
            pdb_file = orig_dir / f"{pdb_id}.pdb"
            if pdb_file.exists():
                try:
                    title = " ".join(l[10:].strip() for l in open(pdb_file) if l.startswith("TITLE"))
                    title = " ".join(title.split())[:120]
                except Exception:
                    title = "N/A"
            else:
                title = "N/A"
            rows.append({"Receptor": rec, "PDB_chain": pdb_chain,
                         "PDB_ID": pdb_id, "Description": title})
    df = pd.DataFrame(rows)
    df = df.drop_duplicates(subset=["Receptor", "PDB_chain"])
    df.to_csv(OUT / "TableS3_offtarget_proteins.csv", index=False)
    print(f"[OK] TableS3: {len(df)} rows")

if __name__ == "__main__":
    tableS1()
    tableS2()
    tableS3()
