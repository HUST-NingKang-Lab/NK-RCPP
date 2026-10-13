#!/usr/bin/env python3
"""
Build the step-1 candidate tables (the PepMLM Colab deliverables) in one
uniform schema, and verify that the candidate numbering matches the numbering
used by every downstream docking input.

Output (step1_peptide_generation/output/):
  NK1R_pepmlm_candidates.csv
  NK2R_pepmlm_candidates.csv
  NK3R_pepmlm_candidates.csv

Columns: peptide_id, sequence, pseudo_perplexity, role
  role = candidate | receptor_matched_control

Caveat recorded in the table itself: the NK2R Colab export was not archived,
so the NK2R pseudo-perplexities are unavailable and only the sequence list was
recovered from the docking input file.

Usage
-----
    python scripts/make_candidate_tables.py --source /path/to/peptide_docking \
                                            --repo   /path/to/NK-RCPP
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

CONTROLS = {
    "NK1R": ("Substance P", "RPKPQQFFGLM"),
    "NK2R": ("Neurokinin A derived", "HKTDSFVGLM"),
    "NK3R": ("Neurokinin B", "DMHDFFVGLM"),
}


def from_colab(csv: Path) -> pd.DataFrame:
    df = pd.read_csv(csv)
    df.columns = [c.strip() for c in df.columns]
    df["sequence"] = df["Binder"].astype(str).str.strip()
    return pd.DataFrame(
        {
            "peptide_id": [f"pep{i}" for i in range(1, len(df) + 1)],
            "sequence": df["sequence"],
            "pseudo_perplexity": df["Pseudo Perplexity"].astype(float),
        }
    )


def from_docking_input(csv: Path) -> pd.DataFrame:
    df = pd.read_csv(csv)
    df["complex_name"] = df["complex_name"].astype(str).str.strip()
    df = df[df["complex_name"] != "pep_known"]
    df["sequence"] = df["peptide_description"].astype(str).str.strip()
    out = pd.DataFrame(
        {
            "peptide_id": df["complex_name"].tolist(),
            "sequence": df["sequence"].tolist(),
        }
    )
    out["pseudo_perplexity"] = pd.NA
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True, help="peptide_docking directory")
    ap.add_argument("--repo", required=True, help="NK-RCPP repository root")
    args = ap.parse_args()

    source = Path(args.source).expanduser().resolve()
    repo = Path(args.repo).expanduser().resolve()
    out = repo / "step1_peptide_generation" / "output"
    out.mkdir(parents=True, exist_ok=True)

    sources = {
        "NK1R": from_colab(source / "nk1r_scripts" / "nk1r.csv"),
        "NK2R": from_docking_input(source / "RAPiDock" / "NK2R_peptide.csv"),
        "NK3R": from_colab(source / "nk3r_scripts" / "nk3r.csv"),
    }

    for receptor, table in sources.items():
        name, seq = CONTROLS[receptor]
        table["role"] = "candidate"
        control = pd.DataFrame(
            [{"peptide_id": "control", "sequence": seq,
              "pseudo_perplexity": pd.NA, "role": "receptor_matched_control"}]
        )
        table = pd.concat([table, control], ignore_index=True)
        table.insert(0, "receptor", receptor)
        table.insert(3, "control_name", name if receptor else "")
        table.loc[table["role"] == "candidate", "control_name"] = ""

        table.to_csv(out / f"{receptor}_pepmlm_candidates.csv", index=False)

        # integrity check against the docking input that drives everything downstream
        # (the control row sits last and is named after the ligand, not "pepN")
        if receptor == "NK2R":
            continue
        dock = pd.read_csv(source / "RAPiDock" / f"{receptor}_peptide.csv")
        dock_seq = dock["peptide_description"].astype(str).str.strip().tolist()[:32]
        colab_seq = table.loc[table["role"] == "candidate", "sequence"].tolist()
        assert colab_seq == dock_seq, f"{receptor}: Colab export and docking input disagree"
        print(f"[ok] {receptor}: {len(colab_seq)} candidates, numbering consistent")


if __name__ == "__main__":
    main()
