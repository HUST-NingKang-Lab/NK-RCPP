#!/usr/bin/env python3
"""
Harvest raw RAPiDock outputs into compact, per-step score tables.

Reads the raw RAPiDock run directory (one folder per complex, each holding
`ref2015_score.csv` with 10 ranked poses) and writes:

  step2_positive_docking/output/positive_scores_long.csv
  step2_positive_docking/output/positive_scores_wide.csv
  step4_negative_docking/output/negative_scores_long.csv
  step6_competition_benchmark/output/competition_positive_scores.csv

Usage
-----
    python scripts/harvest_docking_scores.py --source /path/to/peptide_docking/RAPiDock \
                                             --repo   /path/to/NK-RCPP

For the negative stage the peptide sequence is resolved from the docking-input
CSV rather than from the batch label: on NK2R the batch labels (peptide01..19)
do NOT follow the candidate numbering of the positive stage. Matching on the
label silently mis-assigns off-target scores to the wrong peptide.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

# campaign -> (result dir under results/, control complex name)
POSITIVE_CAMPAIGNS = {
    "NK1R": ("NK1R_peptide_len10", "Substance_P"),
    "NK2R": ("NK2R_peptide_len10", "pep_known"),
    "NK3R": ("NK3R_peptide_len10", "Neurokinin_B"),
    "COMP_TOP10": ("COMP_TOP10_NK2R", None),
    "LENGTH_CONTROL": ("LENGTH_CONTROL", None),
}

# campaign -> (result dir under results/, merged docking-input CSV, per-batch input dir)
# The per-batch dir is the fallback for the sequence: the merged competition CSV
# still carries the pre-fix comp00..comp09 labels, while the batch folders and the
# result directories were renamed to comp01..comp10.
NEGATIVE_CAMPAIGNS = {
    "NK1R": ("NK1R_neg", "NK1R_negative_docking_input.csv", "nk1r_neg_inputs"),
    "NK2R": ("NK2R_neg", "negative_docking_input.csv", "nk2r_neg_inputs"),
    "NK3R": ("NK3R_neg", "NK3R_negative_docking_input.csv", "nk3r_neg_inputs"),
    "COMP_TOP10": ("COMP_TOP10_neg", "COMP_TOP10_negative_docking_input.csv",
                   "comp_top10_neg_inputs"),
}

LABEL_RE = re.compile(r"_(comp\d{2}|peptide\d{2}|known)$")


def read_scores(score_csv: Path) -> pd.DataFrame:
    """Read one ref2015_score.csv -> DataFrame[rank, score]."""
    df = pd.read_csv(score_csv)
    df.columns = [c.strip() for c in df.columns]
    df["rank"] = df["file"].str.extract(r"rank(\d+)").astype(int)
    df["score"] = df["ref2015score"].astype(float)
    return df.sort_values("rank")[["rank", "score"]]


def harvest_positive(source: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for campaign, (subdir, control) in POSITIVE_CAMPAIGNS.items():
        root = source / "results" / subdir
        if not root.is_dir():
            print(f"[skip] missing {root}")
            continue
        for complex_dir in sorted(root.iterdir()):
            score_csv = complex_dir / "ref2015_score.csv"
            if not score_csv.is_file():
                continue
            df = read_scores(score_csv)
            for _, r in df.iterrows():
                rows.append(
                    {
                        "campaign": campaign,
                        "complex_name": complex_dir.name,
                        "is_control": complex_dir.name == control,
                        "rank": int(r["rank"]),
                        "ref2015_score": r["score"],
                    }
                )
    long = pd.DataFrame(rows)
    wide = (
        long.pivot_table(
            index=["campaign", "complex_name", "is_control"],
            columns="rank",
            values="ref2015_score",
        )
        .add_prefix("rank")
        .reset_index()
    )
    wide["best_score"] = wide.filter(like="rank").min(axis=1)
    return long, wide


def harvest_negative(source: Path) -> pd.DataFrame:
    rows = []
    for campaign, (subdir, input_csv, batch_dir) in NEGATIVE_CAMPAIGNS.items():
        root = source / "results" / subdir
        inputs = source / input_csv
        if not root.is_dir() or not inputs.is_file():
            print(f"[skip] missing {root} or {inputs}")
            continue

        meta = pd.read_csv(inputs)
        meta["sequence"] = meta["peptide_description"].astype(str).str.strip()
        meta["target"] = (
            meta["protein_description"].astype(str).str.strip().str.rsplit("/", n=1).str[-1]
            .str.replace(r"\.pdb$", "", regex=True)
        )
        lookup = meta.set_index("complex_name")[["target", "sequence"]].to_dict("index")

        # batch label -> sequence, read straight from the per-batch input files
        by_batch = {}
        bdir = source / batch_dir
        if bdir.is_dir():
            for f in sorted(bdir.glob("*_docking_input.csv")):
                label = f.name.replace("_docking_input.csv", "")
                by_batch[label] = pd.read_csv(f)["peptide_description"].astype(str).str.strip().iloc[0]

        # Two layouts occur in practice: <root>/<batch>/<complex>/ and a flat
        # <root>/<complex>/. Walking for the score file handles both; the batch
        # label is recovered from the complex name (…_peptide07 / …_comp03 / …_known).
        for score_csv in sorted(root.rglob("ref2015_score.csv")):
            complex_dir = score_csv.parent
            m = LABEL_RE.search(complex_dir.name)
            label = m.group(1) if m else None
            # the target is always recoverable from the folder name itself
            target = (complex_dir.name[: -len(label) - 1] if label else complex_dir.name)
            info = lookup.get(complex_dir.name, {})
            # The per-batch input file wins: it is keyed by the batch label that the
            # result folder itself carries, whereas the merged CSV is keyed by complex
            # name and can be stale (its competition rows still use the pre-fix
            # comp00..comp09 labels, which would silently return a neighbouring
            # peptide's sequence).
            seq = by_batch.get(label) or info.get("sequence")
            if seq and info.get("sequence") and seq != info["sequence"]:
                print(f"[warn] {campaign}/{complex_dir.name}: merged CSV says "
                      f"{info['sequence']}, batch input says {seq} -> using batch input")
            df = read_scores(score_csv)
            rows.append(
                {
                    "campaign": campaign,
                    "batch": label or complex_dir.name,
                    "complex_name": complex_dir.name,
                    "target": target,
                    "sequence": seq,
                    "rank1_score": float(df.iloc[0]["score"]),
                    "best_score": float(df["score"].min()),
                    "n_poses": int(len(df)),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True, help="RAPiDock run directory (contains results/)")
    ap.add_argument("--repo", required=True, help="NK-RCPP repository root")
    args = ap.parse_args()

    source = Path(args.source).expanduser().resolve()
    repo = Path(args.repo).expanduser().resolve()

    long, wide = harvest_positive(source)
    out2 = repo / "step2_positive_docking" / "output"
    out2.mkdir(parents=True, exist_ok=True)
    long.to_csv(out2 / "positive_scores_long.csv", index=False)
    wide.to_csv(out2 / "positive_scores_wide.csv", index=False)
    print(f"[ok] positive: {len(long)} pose rows / {long['complex_name'].nunique()} complexes")

    comp = wide[wide["campaign"] == "COMP_TOP10"].copy()
    out6 = repo / "step6_competition_benchmark" / "output"
    out6.mkdir(parents=True, exist_ok=True)
    comp.to_csv(out6 / "competition_positive_scores.csv", index=False)

    neg = harvest_negative(source)
    out4 = repo / "step4_negative_docking" / "output"
    out4.mkdir(parents=True, exist_ok=True)
    neg.to_csv(out4 / "negative_scores_long.csv", index=False)
    if len(neg):
        summary = (
            neg.groupby(["campaign", "batch"])
            .agg(n_targets=("target", "nunique"), mean_score=("rank1_score", "mean"))
            .reset_index()
        )
        summary.to_csv(out4 / "negative_scores_per_batch.csv", index=False)
    print(f"[ok] negative: {len(neg)} complex rows across {neg['batch'].nunique()} batches")


if __name__ == "__main__":
    main()
