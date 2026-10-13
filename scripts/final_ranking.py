#!/usr/bin/env python3
"""
最终综合分析: 合并正向/负向/配对数据, 生成每个受体的最终候选排名
指标:
  ΔPos = pos(candidate) - pos(known)   <0 更好 (结合目标更紧)
  ΔNeg = neg_mean(candidate) - neg_mean(known)  >0 更好 (结合脱靶更弱)
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
# ---------------------------------------------------------------------------


RECEPTORS = {
    "NK1R": {"known_pos": 574.39, "known_neg": 1331.716364, "known_name": "Substance_P"},
    "NK2R": {"known_pos": 771.13, "known_neg": 1049.023704, "known_name": "HKTDSFVGLM"},
    "NK3R": {"known_pos": 2240.54, "known_neg": 1145.922000, "known_name": "Neurokinin_B"},
}

def main():
    for receptor in ["NK1R", "NK2R", "NK3R"]:
        info = RECEPTORS[receptor]
        df = pd.read_csv(ANALYSIS / f"{receptor}_summary.csv")
        pw = pd.read_csv(ANALYSIS / f"{receptor}_pairwise_selectivity.csv")
        # Merge pairwise by label
        # Merge by sequence (robust across numbering schemes)
        merged = df.merge(pw, on="sequence", how="left")

        merged["ΔPos"] = merged["pos_score"] - info["known_pos"]
        merged["ΔNeg"] = merged["neg_mean"] - info["known_neg"]

        # Classification: better binding AND safer
        merged["better_binding"] = merged["ΔPos"] < 0
        merged["safer"] = merged["frac_safer"] > 0.5
        merged["both"] = merged["better_binding"] & merged["safer"]

        merged = merged.sort_values(["both", "ΔPos", "ΔNeg"], ascending=[False, True, False])
        merged.to_csv(ANALYSIS / f"{receptor}_final_ranking.csv", index=False)
        print(f"===== {receptor} FINAL RANKING =====")
        cols = ["peptide_id", "sequence", "pos_score", "ΔPos", "neg_mean", "ΔNeg",
                "frac_safer", "delta_mean", "better_binding", "safer", "both"]
        print(merged[cols].to_string(index=False))
        print()

if __name__ == "__main__":
    main()
