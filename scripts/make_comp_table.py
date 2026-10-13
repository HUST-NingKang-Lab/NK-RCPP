#!/usr/bin/env python3
"""
生成竞赛对比附表 (Table S4) 与差异分析文档数据
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
FIGDIR = BASE / "figures"
TABLES = FIGDIR / "tables"
# ---------------------------------------------------------------------------

TABLES.mkdir(parents=True, exist_ok=True)


def main():
    comp = pd.read_csv(ANALYSIS / "COMP_TOP10_summary.csv")
    ours = pd.read_csv(ANALYSIS / "NK2R_final_ranking.csv")
    ctl = pd.read_csv(ANALYSIS / "length_control.csv")
    stats = ctl.groupby("len")["pos"].agg(["mean", "std"])

    ours["code"] = ["NK2R-K" if p == "pep_known" else f"NK2R-{int(p.replace('pep','')):02d}"
                    for p in ours["peptide_id"]]
    ours["length"] = 10
    ours["z_vs_length"] = (ours["pos_score"] - stats.loc[10, "mean"]) / stats.loc[10, "std"]
    ours["source"] = "This study"
    ours["EC50_nM"] = np.nan
    ours["peptide_name"] = ""

    comp["source"] = "Competition Top10"

    cols = ["source", "code", "peptide_name", "sequence", "length", "EC50_nM",
            "pos_score", "z_vs_length", "ΔPos", "neg_mean", "ΔNeg",
            "frac_safer", "delta_mean", "both"]
    both = pd.concat([comp[[c for c in cols if c in comp.columns]],
                      ours[[c for c in cols if c in ours.columns]]], ignore_index=True)
    both = both.sort_values(["source", "z_vs_length"])
    both.to_csv(TABLES / "TableS4_competition_comparison.csv", index=False)
    print(f"[OK] TableS4: {len(both)} rows")

    # ---- summary statistics for the manuscript ----
    print("\n=== 关键统计 ===")
    print(f"竞赛Top10: z<0 (优于同长随机对照): {(comp['z_vs_length']<0).sum()}/10")
    print(f"竞赛Top10: ΔPos<0 (优于已知肽):   {(comp['ΔPos']<0).sum()}/10")
    print(f"竞赛Top10: frac_safer>0.5:          {(comp['frac_safer']>0.5).sum()}/10")
    print(f"竞赛Top10: 双优(both):              {int(comp['both'].sum())}/10")
    print()
    oo = ours.dropna(subset=["frac_safer"])
    print(f"本研究(18个有负向数据): 双优:       {int(oo['both'].sum())}/{len(oo)}")
    print(f"本研究(32个全部): ΔPos<0:           {(ours['ΔPos']<0).sum()}/32")
    print(f"本研究(32个全部): z_vs_length<0:    {(ours['z_vs_length']<0).sum()}/32")

    # length-matched subgroup (9-mers from competition vs our 10-mers)
    c9 = comp[comp["length"] == 9]
    print(f"\n9-mer 竞赛肽 (长度最接近): 双优 {int(c9['both'].sum())}/{len(c9)}")
    print(f"  {c9[['code','sequence','ΔPos','frac_safer','both']].to_string(index=False)}")

    # correlation
    from scipy import stats as st
    s2 = comp.dropna(subset=["EC50_nM", "z_vs_length"])
    r = st.spearmanr(s2["EC50_nM"], s2["z_vs_length"])
    print(f"\nEC50 vs z_vs_length: ρ={r.statistic:+.2f}, P={r.pvalue:.3f}, n={len(s2)}")
    # rank comparison of EC50 vs docking rank (top3 overlap)
    comp_sorted_ec50 = comp.dropna(subset=["EC50_nM"]).sort_values("EC50_nM")["code"].tolist()
    comp_sorted_z = comp.dropna(subset=["EC50_nM"]).sort_values("z_vs_length")["code"].tolist()
    print(f"EC50 排名前3:    {comp_sorted_ec50[:3]}")
    print(f"对接z 排名前3:   {comp_sorted_z[:3]}")


if __name__ == "__main__":
    main()
