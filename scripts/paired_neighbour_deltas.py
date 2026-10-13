#!/usr/bin/env python3
"""
Rebuild the pair-level matched structural-neighbour comparisons from
step 4's per-complex scores.

The manuscript's audit had to treat the per-candidate aggregates as *inherited*
results, because the packaged archive it reanalysed contained only the
aggregates (no pair-level rows). This repository keeps the per-complex REF15
scores, so the comparison behind every aggregate can be recomputed here:

    for each structural neighbour s of a receptor panel
        delta(s) = candidate_score(s) - receptor_matched_control_score(s)
        n_higher = #{ s : delta(s) > 0 }
        fraction = n_higher / n_common

That reproduces the `n_common` / `n_safer_delta_pos` / `frac_safer` columns of
step 5's `*_pairwise_selectivity.csv` from the raw numbers.

Outputs
-------
step5_analysis/output/paired_neighbour_deltas.csv   one row per (receptor, candidate, neighbour)
step5_analysis/output/paired_deltas_summary.csv     one row per (receptor, candidate), with the
                                                    recomputed counts next to the shipped ones

Usage
-----
    python scripts/paired_neighbour_deltas.py --repo .
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

RECEPTORS = ["NK1R", "NK2R", "NK3R"]

# campaign in negative_scores_long.csv -> receptor, and the batch that holds the control
CONTROL_BATCH = {"NK1R": "known", "NK2R": "peptide19", "NK3R": "known"}
CAMPAIGN = {"NK1R": "NK1R", "NK2R": "NK2R", "NK3R": "NK3R"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=".")
    args = ap.parse_args()
    repo = Path(args.repo).expanduser().resolve()

    neg = pd.read_csv(repo / "step4_negative_docking" / "output" / "negative_scores_long.csv")
    out = repo / "step5_analysis" / "output"

    rows, summaries = [], []
    for receptor in RECEPTORS:
        sub = neg[neg["campaign"] == CAMPAIGN[receptor]]
        control_seq = sub[sub["batch"] == CONTROL_BATCH[receptor]]["sequence"].dropna()
        if control_seq.empty:
            print(f"[skip] {receptor}: no control batch found")
            continue
        control_seq = control_seq.iloc[0]

        # candidate batch -> sequence; the control's own sequence must be dropped
        batches = sub.groupby("batch")["sequence"].first()
        cand_batches = batches[(batches != control_seq) & batches.notna()]

        ctrl = (
            sub[(sub["batch"] == CONTROL_BATCH[receptor])]
            .set_index("target")["rank1_score"]
        )

        wide = (
            sub[sub["batch"].isin(cand_batches.index)]
            .pivot_table(index="target", columns="batch", values="rank1_score")
        )
        for batch in wide.columns:
            seq = cand_batches[batch]
            col = wide[batch].dropna()
            common = col.index.intersection(ctrl.index)
            if not len(common):
                continue
            delta = col.loc[common] - ctrl.loc[common]
            n_common = int(len(delta))
            n_higher = int((delta > 0).sum())

            for target, d in delta.items():
                rows.append(
                    {
                        "receptor": receptor,
                        "batch": batch,
                        "sequence": seq,
                        "neighbour": target,
                        "candidate_score": float(col.loc[target]),
                        "control_score": float(ctrl.loc[target]),
                        "delta": float(d),
                        "higher_than_control": bool(d > 0),
                    }
                )
            summaries.append(
                {
                    "receptor": receptor,
                    "batch": batch,
                    "sequence": seq,
                    "n_common": n_common,
                    "n_higher": n_higher,
                    "fraction": n_higher / n_common,
                    "delta_mean": float(delta.mean()),
                    "delta_median": float(delta.median()),
                }
            )

    pairs = pd.DataFrame(rows)
    summ = pd.DataFrame(summaries)
    pairs.to_csv(out / "paired_neighbour_deltas.csv", index=False)
    summ.to_csv(out / "paired_deltas_summary.csv", index=False)
    print(f"[ok] {len(pairs)} pair-level comparisons over {len(summ)} candidates")

    # cross-check against the shipped aggregates, by sequence (never by batch label)
    for receptor in RECEPTORS:
        ship_path = out / f"{receptor}_pairwise_selectivity.csv"
        if not ship_path.is_file():
            continue
        ship = pd.read_csv(ship_path)
        mine = summ[summ["receptor"] == receptor]
        merged = ship.merge(mine, on="sequence", suffixes=("_shipped", "_recomputed"))
        if merged.empty:
            continue
        agree = (merged["n_common_shipped"] == merged["n_common_recomputed"]) & (
            merged["n_safer_delta_pos"] == merged["n_higher"]
        )
        n_bad = int((~agree).sum())
        print(f"[check] {receptor}: {len(merged)} candidates compared, "
              f"{n_bad} disagree with the shipped aggregates")
        if n_bad:
            print(merged.loc[~agree, ["sequence", "n_common_shipped",
                                      "n_common_recomputed", "n_safer_delta_pos",
                                      "n_higher"]].to_string())


if __name__ == "__main__":
    main()
