#!/usr/bin/env python3
"""
Condense the Foldseek 3Di alignments into a per-receptor hit table.

`alis_pdb100.m8` from the Foldseek web API carries 21 columns per alignment
(the 12 BLAST columns plus probability/e-value/bitscore and the aligned
sequences). Only the best alignment per structural neighbour is kept, so the
table stays small enough to ship while remaining a faithful summary of what
the search returned.

Writes step3_foldseek_search/output/<RECEPTOR>/pdb100_best_hits.csv and copies
the small companion outputs (taxonomy report, hit PDB id list, CD-HIT-50
representative FASTA) next to it.

Usage
-----
    python scripts/harvest_foldseek_hits.py --source /path/to/peptide_docking \
                                            --repo   /path/to/NK-RCPP
"""
from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

import pandas as pd

RECEPTORS = {"NK1R": "nk1r_foldseek", "NK2R": "nk2r_foldseek", "NK3R": "nk3r_foldseek"}

# the NK2R panel was generated before the per-receptor file naming was adopted
PANEL_CSV = {
    "NK1R": "NK1R_negative_docking_input.csv",
    "NK2R": "negative_docking_input.csv",
    "NK3R": "NK3R_negative_docking_input.csv",
}

# 21-column placement delivered by the Foldseek API for 3diaa/pdb100 searches
COLS = [
    "query", "target", "fident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "prob", "evalue", "bits",
    "c14", "c15", "qseq", "tseq", "tcoords", "tseq2", "taxid", "taxname",
]

TARGET_RE = re.compile(r"^(?P<pdb>[0-9a-z]{4})(?:-assembly\d+)?\.cif\.gz_(?P<chain>\S+)$")


def parse_target(raw: str) -> tuple[str, str]:
    """'6e59-assembly1.cif.gz_A' -> ('6e59', 'A').

    A repeated chain is suffixed by foldseek ('…_A-2' = second alignment of
    chain A); the suffix is stripped so the key matches the panel folder name.
    """
    m = TARGET_RE.match(raw.split()[0])
    if not m:
        return raw.split()[0], ""
    return m.group("pdb"), re.sub(r"-\d+$", "", m.group("chain"))


def harvest_one(fs_dir: Path, panel: set[str]) -> pd.DataFrame:
    m8 = fs_dir / "alis_pdb100.m8"
    df = pd.read_csv(m8, sep="\t", header=None, names=COLS, low_memory=False)
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
    df["bits"] = pd.to_numeric(df["bits"], errors="coerce")

    parsed = df["target"].map(parse_target)
    df["pdb"] = [p[0] for p in parsed]
    df["chain"] = [p[1] for p in parsed]
    df["neighbour"] = df["pdb"] + "_" + df["chain"]
    df["description"] = df["target"].str.replace(r"^\S+\s*", "", regex=True)

    # best alignment per structural neighbour (largest bitscore)
    best = (
        df.sort_values("bits", ascending=False)
        .drop_duplicates(subset="neighbour", keep="first")
        .copy()
    )
    best["in_negative_panel"] = best["neighbour"].isin(panel)
    keep = ["query", "neighbour", "pdb", "chain", "description", "fident",
            "alnlen", "evalue", "bits", "taxname", "in_negative_panel"]
    return (
        best[keep]
        .sort_values("bits", ascending=False)
        .rename(columns={"fident": "seq_identity_pct", "alnlen": "aln_len"})
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True, help="peptide_docking directory")
    ap.add_argument("--repo", required=True, help="NK-RCPP repository root")
    args = ap.parse_args()

    source = Path(args.source).expanduser().resolve()
    repo = Path(args.repo).expanduser().resolve()
    out_root = repo / "step3_foldseek_search" / "output"

    for receptor, subdir in RECEPTORS.items():
        fs_dir = source / subdir
        if not fs_dir.is_dir():
            print(f"[skip] {fs_dir}")
            continue
        out = out_root / receptor
        out.mkdir(parents=True, exist_ok=True)

        # structural neighbours that survived into the negative-docking panel
        panel_csv = source / "RAPiDock" / PANEL_CSV[receptor]
        panel = set()
        if panel_csv.is_file():
            p = pd.read_csv(panel_csv)
            panel = set(
                p["protein_description"].astype(str).apply(lambda s: Path(s).stem)
            )

        hits = harvest_one(fs_dir, panel)
        hits.to_csv(out / "pdb100_best_hits.csv", index=False)

        missing = panel - set(hits["neighbour"])
        if missing:
            print(f"[note] {receptor}: {len(missing)} panel entries absent from the "
                  f"pdb100 hit list: {sorted(missing)[:5]}")

        for name in ["alis_pdb100_report.m8", "pdb_ids.txt", "alis_pdb100_cdhit50.fasta"]:
            src = fs_dir / name
            if src.is_file():
                shutil.copy2(src, out / name)

        n_panel = int(hits["in_negative_panel"].sum())
        print(f"[ok] {receptor}: {len(hits)} neighbours, {n_panel} used in negative panel")


if __name__ == "__main__":
    main()
