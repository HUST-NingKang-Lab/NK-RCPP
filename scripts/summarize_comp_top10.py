#!/usr/bin/env python3
"""
汇总竞赛 Top10 肽的对接结果，并与我们的 NK2R 候选对比。
关键: 加入长度对照校正 (随机组成匹配肽的基线)，因 ref2015 分数受肽长强烈影响。

输出:
  analysis/COMP_TOP10_summary.csv       竞赛肽结果
  analysis/NK2R_comparison.csv          合并对比表
  analysis/length_control.csv           长度对照原始数据
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

RESULTS = RESULTS_DIR
KNOWN_POS = 771.13
KNOWN_SEQ = "HKTDSFVGLM"


def read_best(score_csv):
    try:
        return float(pd.read_csv(score_csv)["ref2015score"].min())
    except Exception:
        return np.nan


def load_length_control():
    """Random composition-matched peptides docked to NK2R, by length."""
    base = RESULTS / "LENGTH_CONTROL"
    rows = []
    for d in sorted(base.glob("ctl*mer_*")):
        name = d.name
        L = int(name.split("mer")[0].replace("ctl", ""))
        rows.append({"name": name, "len": L, "pos": read_best(d / "ref2015_score.csv")})
    return pd.DataFrame(rows).dropna()


def length_stats(ctl):
    """Per-length mean/std of random control scores."""
    return ctl.groupby("len")["pos"].agg(["mean", "std", "count"])


def zscore_by_length(score, length, stats):
    """Z-score of a peptide score vs length-matched random controls."""
    if length not in stats.index:
        return np.nan
    m, s = stats.loc[length, "mean"], stats.loc[length, "std"]
    return (score - m) / s if s > 0 else np.nan


def main():
    ctl = load_length_control()
    ctl.to_csv(ANALYSIS / "length_control.csv", index=False)
    stats = length_stats(ctl)
    print("=== 长度对照 (随机组成匹配肽 vs NK2R) ===")
    print(stats.round(1).to_string())
    print()

    meta = pd.read_csv(ANALYSIS / "competition_top10.csv")
    NEG_ROOT = RESULTS / "COMP_TOP10_neg"

    # ---- 1. Competition Top10 ----
    rows = []
    for i, r in meta.iterrows():
        pid = f"comp{i+1:02d}"
        pos = read_best(RESULTS / "COMP_TOP10_NK2R" / pid / "ref2015_score.csv")
        # negative results live in dirs named <protein>_<pid>, e.g. 3oe6_A_comp01
        neg = {}
        for d in NEG_ROOT.glob(f"*_{pid}"):
            sp = d / "ref2015_score.csv"
            if sp.exists():
                neg[d.name[: -len(pid) - 1]] = read_best(sp)
        rows.append({
            "code": f"COMP-{i+1:02d}", "peptide_name": r["peptide_name"],
            "sequence": r["sequence"], "EC50_nM": r["EC50_nM"], "rank": r["rank"],
            "length": len(str(r["sequence"])), "pos_score": pos,
            "n_neg": len(neg), "_neg": neg,
        })

    # known ligand negative scores (own batch for paired comparison)
    known_neg = {}
    for d in NEG_ROOT.glob("*_known"):
        sp = d / "ref2015_score.csv"
        if sp.exists():
            known_neg[d.name[:-len("_known")]] = read_best(sp)
    known_neg_mean = np.mean(list(known_neg.values())) if known_neg else np.nan

    comp_rows = []
    for row in rows:
        d = row.pop("_neg")
        common = set(d) & set(known_neg)
        deltas = {p: d[p] - known_neg[p] for p in common} if common else {}
        n_safer = sum(1 for v in deltas.values() if v > 0)
        vals = [d[p] for p in d if not np.isnan(d[p])]
        row["neg_mean"] = np.mean(vals) if vals else np.nan
        row["ΔPos"] = row["pos_score"] - KNOWN_POS if pd.notna(row["pos_score"]) else np.nan
        row["ΔNeg"] = (row["neg_mean"] - known_neg_mean) if row["neg_mean"] == row["neg_mean"] else np.nan
        row["frac_safer"] = (n_safer / len(common)) if common else np.nan
        row["delta_mean"] = np.mean(list(deltas.values())) if deltas else np.nan
        row["both"] = ((row["ΔPos"] < 0) and (row["frac_safer"] > 0.5)) if pd.notna(row["frac_safer"]) else np.nan
        # length-corrected
        row["z_vs_length"] = zscore_by_length(row["pos_score"], row["length"], stats)
        comp_rows.append(row)

    comp = pd.DataFrame(comp_rows).sort_values("z_vs_length")
    comp.to_csv(ANALYSIS / "COMP_TOP10_summary.csv", index=False)

    print("=== 竞赛 Top10 对接结果 (按长度校正 z 值排序) ===")
    print(comp[["code", "peptide_name", "sequence", "length", "EC50_nM", "pos_score",
                "z_vs_length", "ΔPos", "frac_safer", "both"]].round(3).to_string(index=False))
    print(f"\n[known] pos={KNOWN_POS}  z_vs_length={zscore_by_length(KNOWN_POS,10,stats):+.2f}  "
          f"neg_mean={known_neg_mean:.1f} (n={len(known_neg)})")

    # ---- 2. Our candidates + length-corrected ----
    ours = pd.read_csv(ANALYSIS / "NK2R_final_ranking.csv")
    ours["code"] = ["NK2R-K" if p == "pep_known" else f"NK2R-{int(p.replace('pep','')):02d}"
                    for p in ours["peptide_id"]]
    ours["length"] = 10
    ours["z_vs_length"] = ours["pos_score"].apply(lambda s: zscore_by_length(s, 10, stats))

    # ---- 3. Combined ----
    cols = ["code", "sequence", "length", "pos_score", "z_vs_length", "ΔPos",
            "neg_mean", "ΔNeg", "frac_safer", "delta_mean", "both"]
    comp_out = comp[cols].copy(); comp_out["source"] = "Competition Top10"
    comp_out["EC50_nM"] = comp["EC50_nM"].values
    ours_out = ours[cols].copy(); ours_out["source"] = "This study"
    ours_out["EC50_nM"] = np.nan
    combined = pd.concat([comp_out, ours_out], ignore_index=True).sort_values("z_vs_length")
    combined.to_csv(ANALYSIS / "NK2R_comparison.csv", index=False)
    print(f"\n[OK] {ANALYSIS}/COMP_TOP10_summary.csv, NK2R_comparison.csv")

    # ---- 4. Stats ----
    print("\n=== 对比摘要 ===")
    print(f"竞赛Top10 中 z_vs_length < 0 (优于同长随机): {int((comp['z_vs_length']<0).sum())}/10")
    print(f"竞赛Top10 中 ΔPos < 0 (优于已知肽): {int((comp['ΔPos']<0).sum())}/10")
    print(f"竞赛Top10 中双优: {int(comp['both'].sum())}/10")
    oc = ours.dropna(subset=["frac_safer"])
    print(f"我方候选双优: {int(oc['both'].sum())}/{len(oc)} (有负向数据的)")
    print(f"我方 z_vs_length < 0: {int((ours['z_vs_length']<0).sum())}/32")

    # ---- 5. EC50 correlation with docking ----
    from scipy import stats as st
    s2 = comp.dropna(subset=["EC50_nM", "pos_score"])
    r1 = st.spearmanr(s2["EC50_nM"], s2["pos_score"])
    r2 = st.spearmanr(s2["EC50_nM"], s2["z_vs_length"])
    print(f"\nEC50 vs 原始 pos_score:  ρ={r1.statistic:+.2f} (P={r1.pvalue:.3f})")
    print(f"EC50 vs 长度校正 z 值:   ρ={r2.statistic:+.2f} (P={r2.pvalue:.3f})")


if __name__ == "__main__":
    main()
