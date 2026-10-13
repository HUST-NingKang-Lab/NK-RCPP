#!/usr/bin/env python3
"""
Per-protein 配对选择性分析
对每个脱靶蛋白: delta = 候选肽分数 - 已知肽分数 (同一蛋白上)
delta > 0 表示候选肽对该蛋白结合比已知肽更弱 (更安全)
同时比较正向分数: 候选 vs 已知 (越低越好)
"""
import pandas as pd
import os
import numpy as np
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
RESULTS_DIR = RAPIDOCK / "results"
# ---------------------------------------------------------------------------

OUTPUT_DIR = ANALYSIS
OUTPUT_DIR.mkdir(exist_ok=True)

RECEPTORS = {
    "NK1R": {"neg_dir": "NK1R_neg", "known": "Substance_P", "known_seq": "RPKPQQFFGLM"},
    "NK2R": {"neg_dir": "NK2R_neg", "known": "pep_known", "known_seq": "HKTDSFVGLM"},
    "NK3R": {"neg_dir": "NK3R_neg", "known": "Neurokinin_B", "known_seq": "DMHDFFVGLM"},
}

def read_best_score(score_csv):
    try:
        df = pd.read_csv(score_csv)
        return float(df["ref2015score"].min())
    except Exception:
        return np.nan

def resolve_label_seq(receptor, label):
    """Resolve neg label (peptideNN/known) -> sequence from neg input CSVs."""
    neg_inputs = RESULTS_DIR.parent / f"{receptor.lower()}_neg_inputs"
    if label == "known":
        f = neg_inputs / "known_peptide_docking_input.csv"
        if f.exists():
            return str(pd.read_csv(f)["peptide_description"].iloc[0]).strip()
        return None
    f = neg_inputs / f"peptide{label}_docking_input.csv"
    if f.exists():
        return str(pd.read_csv(f)["peptide_description"].iloc[0]).strip()
    return None

def build_protein_scores(receptor):
    """Return dict: peptide_seq -> {protein_id: best_score}"""
    info = RECEPTORS[receptor]
    neg_dir = RESULTS_DIR / info["neg_dir"]
    seq_scores = {}
    for score_file in neg_dir.glob("*/**/ref2015_score.csv"):
        complex_name = score_file.parent.name
        if "_peptide" in complex_name:
            label = complex_name.split("_peptide")[-1]
        elif "_known" in complex_name:
            label = "known"
        else:
            continue
        protein = complex_name.rsplit("_", 1)[0]
        score = read_best_score(score_file)
        seq_scores.setdefault(label, {})[protein] = score
    # Resolve labels to sequences, merge by sequence
    by_seq = {}
    for label, scores in seq_scores.items():
        seq = resolve_label_seq(receptor, label)
        if seq:
            by_seq.setdefault(seq, {}).update(scores)
    return by_seq

def main():
    for receptor in ["NK1R", "NK2R", "NK3R"]:
        info = RECEPTORS[receptor]
        seq_scores = build_protein_scores(receptor)
        # known is now keyed by sequence
        known_seq = info["known_seq"]
        if known_seq not in seq_scores:
            print(f"[WARN] {receptor}: known sequence not found. Sequences: {list(seq_scores.keys())[:5]}...")
            continue
        known_scores = seq_scores.get(known_seq, {})
        print(f"=== {receptor}: known peptide ({known_seq}) neg scores on {len(known_scores)} proteins ===")

        rows = []
        for seq, scores in seq_scores.items():
            if seq == known_seq:
                continue
            common = set(scores) & set(known_scores)
            if not common:
                continue
            # For each protein, delta = candidate - known (positive = candidate binds weaker = safer)
            deltas = {p: scores[p] - known_scores[p] for p in common}
            n_positive = sum(1 for d in deltas.values() if d > 0)
            rows.append({
                "sequence": seq,
                "n_common": len(common),
                "n_safer_delta_pos": n_positive,
                "frac_safer": n_positive / len(common),
                "delta_mean": np.mean(list(deltas.values())),
                "delta_median": np.median(list(deltas.values())),
                "delta_min": min(deltas.values()),
            })
        df = pd.DataFrame(rows).sort_values("delta_mean", ascending=False)
        df.to_csv(OUTPUT_DIR / f"{receptor}_pairwise_selectivity.csv", index=False)
        print(f"  {len(df)} peptides analyzed")
        print(df.head(5).to_string(index=False))
        print()

if __name__ == "__main__":
    main()
