#!/usr/bin/env python3
"""
长度对照实验: 生成组成匹配的随机肽 (6/9/10/15-mer)，对接 NK2R
用于量化 ref2015 分数的长度偏差，使竞赛肽与本研究的 10-mer 可比。
每长度 5 条随机肽，氨基酸组成与竞赛 Top10 对应长度组的组成匹配。
"""
import random
import os
import pandas as pd
from pathlib import Path

# --- paths ---------------------------------------------------------------
# Repository-relative: the scripts run from anywhere inside a checkout.
# Raw RAPiDock run directories and the competition workbook are not shipped;
# set NK_RCPP_RAW / NK_RCPP_COMP_XLSX (or NK_RCPP_BASE) to point at them.
# See scripts/README.md.
BASE = Path(os.environ.get("NK_RCPP_BASE", Path(__file__).resolve().parent.parent))
RAW = Path(os.environ.get("NK_RCPP_RAW", BASE / "raw"))
DOCKING = RAW / "peptide_docking"
RAPIDOCK = DOCKING / "RAPiDock"
# ---------------------------------------------------------------------------

random.seed(42)
NK2R_PDB = "peptide_docking/pdbs/NK2R.pdb"
OUT_CSV = RAPIDOCK / "LENGTH_CONTROL_NK2R_peptide.csv"
OUT_DIR = RAPIDOCK / "length_control_inputs"

# Composition-matched pools derived from competition Top10 by length group
POOLS = {
    6:  list("CFYFLM"),                                    # from COMP-03
    9:  list("GFPCFYFLMYCMEGAKMMPFPCFYFLMRFPCFYFLMAFPCFYFLM"),  # 5× 9-mers
    10: list("RSPRQCCFQMHKTDSFVGLM"),                      # our NK2R candidates ref
    15: list("CPYTKLNRCCFYFLMCQRFRPGKCCFYFLMCQTPAPGKCCFYFLMCMTPQPGKCCFYFLM"),
}
N_PER_LENGTH = 5


def main():
    rows = []
    for L, pool in POOLS.items():
        for i in range(N_PER_LENGTH):
            seq = "".join(random.choice(pool) for _ in range(L))
            rows.append({"complex_name": f"ctl{L}mer_{i+1:02d}",
                         "protein_description": NK2R_PDB,
                         "peptide_description": seq})
    df = pd.DataFrame(rows, columns=["complex_name", "protein_description", "peptide_description"])
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for _, r in df.iterrows():
        sub = df[df["complex_name"] == r["complex_name"]]
        sub.to_csv(OUT_DIR / f'{r["complex_name"]}_docking_input.csv', index=False)
    print(f"[OK] {len(df)} control peptides -> {OUT_CSV}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
