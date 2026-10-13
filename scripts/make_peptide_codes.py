#!/usr/bin/env python3
"""
生成论文肽段代号 → 序列映射表 (Supplementary Table)
代号格式: {NK1R,NK2R,NK3R}-01..32, 已知肽 = {NK1R,NK2R,NK3R}-K
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

KNOWN = {"NK1R": ("Substance_P", "RPKPQQFFGLM"),
         "NK2R": ("pep_known", "HKTDSFVGLM"),
         "NK3R": ("Neurokinin_B", "DMHDFFVGLM")}

def main():
    rows = []
    for rec in ["NK1R", "NK2R", "NK3R"]:
        df = pd.read_csv(ANALYSIS / f"{rec}_summary.csv")
        for _, r in df.iterrows():
            pid = r["peptide_id"]
            seq = str(r["sequence"]).strip()
            if pid in ["Substance_P", "pep_known", "Neurokinin_B"]:
                code = f"{rec}-K"
            else:
                num = pid.replace("pep", "")
                code = f"{rec}-{int(num):02d}"
            rows.append({
                "Code": code,
                "Receptor": rec,
                "PepMLM_ID": pid,
                "Sequence": seq,
                "Pseudo_Perplexity": r.get("pseudo_pp", np.nan),
            })
    df = pd.DataFrame(rows)
    df.to_csv(ANALYSIS / "peptide_codes.csv", index=False)
    print(df.to_string(index=False))

if __name__ == "__main__":
    main()
