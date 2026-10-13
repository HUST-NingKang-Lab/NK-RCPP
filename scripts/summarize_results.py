#!/usr/bin/env python3
"""
汇总三个受体的正向/负向对接结果，计算选择性指标
输出: analysis/ 下的 CSV 文件

注意: 负向数据按"序列"匹配，因为 NK2R 的 peptideNN 编号与 pepN 编号不对应。
"""
import os
import glob
import pandas as pd
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
    "NK1R": {"pos_dir": "NK1R_peptide_len10", "neg_dir": "NK1R_neg",
             "known": "Substance_P", "known_seq": "RPKPQQFFGLM"},
    "NK2R": {"pos_dir": "NK2R_peptide_len10", "neg_dir": "NK2R_neg",
             "known": "pep_known", "known_seq": "HKTDSFVGLM"},
    "NK3R": {"pos_dir": "NK3R_peptide_len10", "neg_dir": "NK3R_neg",
             "known": "Neurokinin_B", "known_seq": "DMHDFFVGLM"},
}


def read_best_score(score_csv: Path) -> float:
    """Read ref2015_score.csv, return best (lowest) score."""
    try:
        df = pd.read_csv(score_csv)
        return float(df["ref2015score"].min())
    except Exception:
        return np.nan


def read_peptide_map(receptor):
    """Read the positive docking CSV to map complex_name -> peptide sequence."""
    pos_csv = RESULTS_DIR.parent / f"{receptor}_peptide.csv"
    df = pd.read_csv(pos_csv)
    return {str(k).strip(): str(v).strip() for k, v in zip(df["complex_name"], df["peptide_description"])}


def build_neg_index(receptor):
    """
    Build negative docking index: peptide_label -> list of per-protein best scores.
    peptide_label comes from complex names like 3oe6_A_peptide01 -> peptide01, xxx_known -> known.
    """
    info = RECEPTORS[receptor]
    neg_dir = RESULTS_DIR / info["neg_dir"]
    neg_index = {}
    if not neg_dir.exists():
        return neg_index
    for score_file in neg_dir.glob("*/**/ref2015_score.csv"):
        complex_name = score_file.parent.name
        if "_peptide" in complex_name:
            label = "peptide" + complex_name.split("_peptide")[-1]
        elif "_known" in complex_name:
            label = "known"
        else:
            continue
        score = read_best_score(score_file)
        neg_index.setdefault(label, []).append(score)
    # deduplicate per complex: same complex may appear once; keep as is (scores per protein)
    return neg_index


def resolve_neg_seq(receptor, neg_label):
    """
    Resolve neg label (peptideNN/known) -> sequence using the negative input CSVs.
    Returns sequence or None.
    """
    info = RECEPTORS[receptor]
    neg_inputs_dir = RESULTS_DIR.parent / f"{receptor.lower()}_neg_inputs"
    known_map = {
        "NK1R": "known_peptide_docking_input.csv",
        "NK2R": None,  # NK2R known was peptide19
        "NK3R": "known_peptide_docking_input.csv",
    }
    if neg_label == "known" and known_map.get(receptor):
        f = neg_inputs_dir / known_map[receptor]
        if f.exists():
            df = pd.read_csv(f)
            return str(df["peptide_description"].iloc[0]).strip()
    # peptideNN
    f = neg_inputs_dir / f"{neg_label}_docking_input.csv"
    if f.exists():
        df = pd.read_csv(f)
        return str(df["peptide_description"].iloc[0]).strip()
    # NK2R: also try nk2r_neg_inputs in RAPiDock dir
    f = RESULTS_DIR.parent / f"{receptor.lower()}_neg_inputs" / f"{neg_label}_docking_input.csv"
    if f.exists():
        df = pd.read_csv(f)
        return str(df["peptide_description"].iloc[0]).strip()
    return None


def summarize_receptor(receptor):
    info = RECEPTORS[receptor]
    pos_dir = RESULTS_DIR / info["pos_dir"]
    neg_dir = RESULTS_DIR / info["neg_dir"]

    peptide_map = read_peptide_map(receptor)
    neg_index = build_neg_index(receptor)

    # Resolve each neg label to a sequence
    neg_seq_map = {}  # seq -> list of neg scores
    for label, scores in neg_index.items():
        seq = resolve_neg_seq(receptor, label)
        if seq:
            neg_seq_map.setdefault(seq, []).extend(scores)

    rows = []
    seen = set()

    for pep_dir in sorted(pos_dir.glob("pep*")):
        pep_name = pep_dir.name
        if pep_name in seen:
            continue
        seen.add(pep_name)
        score_csv = pep_dir / "ref2015_score.csv"
        pos_score = read_best_score(score_csv)
        seq = peptide_map.get(pep_name, "")

        neg_scores = neg_seq_map.get(seq, [])
        neg_scores = [s for s in neg_scores if not np.isnan(s)]

        rows.append({
            "receptor": receptor,
            "peptide_id": pep_name,
            "sequence": seq,
            "pos_score": pos_score,
            "neg_mean": np.mean(neg_scores) if neg_scores else np.nan,
            "neg_median": np.median(neg_scores) if neg_scores else np.nan,
            "neg_min": np.min(neg_scores) if neg_scores else np.nan,
            "neg_max": np.max(neg_scores) if neg_scores else np.nan,
            "n_neg": len(neg_scores),
            "selectivity": (np.mean(neg_scores) - pos_score) if neg_scores else np.nan,
        })

    # Add known peptide if not already present
    known_name = info["known"]
    if known_name not in seen:
        score_csv = pos_dir / known_name / "ref2015_score.csv"
        known_pos = read_best_score(score_csv)
        seq = info["known_seq"]
        neg_scores = neg_seq_map.get(seq, [])
        neg_scores = [s for s in neg_scores if not np.isnan(s)]
        rows.append({
            "receptor": receptor,
            "peptide_id": known_name,
            "sequence": seq,
            "pos_score": known_pos,
            "neg_mean": np.mean(neg_scores) if neg_scores else np.nan,
            "neg_median": np.median(neg_scores) if neg_scores else np.nan,
            "neg_min": np.min(neg_scores) if neg_scores else np.nan,
            "neg_max": np.max(neg_scores) if neg_scores else np.nan,
            "n_neg": len(neg_scores),
            "selectivity": (np.mean(neg_scores) - known_pos) if neg_scores else np.nan,
        })

    return pd.DataFrame(rows)


def main():
    all_dfs = []
    for receptor in ["NK1R", "NK2R", "NK3R"]:
        df = summarize_receptor(receptor)
        all_dfs.append(df)
        df.to_csv(OUTPUT_DIR / f"{receptor}_summary.csv", index=False)
        print(f"[{receptor}] {len(df)} peptides summarized")
        n_missing = (df["n_neg"] == 0).sum()
        print(f"  - {n_missing} peptides missing negative data")
        print(df[["peptide_id", "sequence", "pos_score", "neg_mean", "n_neg", "selectivity"]].head(3).to_string())
        print("  ...")

    combined = pd.concat(all_dfs, ignore_index=True)
    combined.to_csv(OUTPUT_DIR / "all_receptors_summary.csv", index=False)
    print(f"\n[ALL] Combined {len(combined)} rows -> {OUTPUT_DIR / 'all_receptors_summary.csv'}")


if __name__ == "__main__":
    main()
