#!/usr/bin/env python3
"""
UMAP 数据准备: 竞赛排行榜前1000 + 本研究全部候选肽
输出: analysis/umap_input.csv  (seq, group, label, length, ...)
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
XLSX = Path(os.environ.get("NK_RCPP_COMP_XLSX",
                  BASE / "step6_competition_benchmark" / "input" /
                  "2025_pep_design_comp_all_entries_with_synthesize_highlight.xlsx"))
# ---------------------------------------------------------------------------

ANALYSIS.mkdir(exist_ok=True)
TOP_N = 1000

# 非标准氨基酸 -> X (ESM 处理为 unknown)
NONSTD = set("OBZUX")


def clean(seq):
    s = str(seq).strip().upper().replace(" ", "").replace("\t", "")
    return "".join("X" if c in NONSTD else c for c in s)


def main():
    rows = []

    # ---- 1. 竞赛前 1000 ----
    df = pd.read_excel(XLSX, sheet_name="2025_pep_design_comp_all_entrie")
    top = df.head(TOP_N).copy()
    for i, r in top.iterrows():
        seq = clean(r["Amino acid Sequence"])
        team = str(r["Team Name"])
        rows.append({
            "seq": seq,
            "group": "leaderboard",
            "label": f"LB{i+1:04d}",
            "rank": i + 1,
            "length": len(seq),
            "ec50": r["EC50 Value (nM)"],
            "name": str(r["Peptide Name"]),
            "is_top10": i < 10,
            "is_own_team": ("zhouyuan" in team.lower()),
            "team": team,
        })
    print(f"[OK] leaderboard top{TOP_N}: {len(top)} 条 "
          f"(其中 zhouyuan_team {sum(1 for r in rows if r.get('is_own_team'))} 条)")

    # ---- 2. 本研究候选 (三个受体) ----
    lb_seqs = {r["seq"] for r in rows}
    for rec in ["NK2R", "NK1R", "NK3R"]:
        summ = pd.read_csv(ANALYSIS / f"{rec}_summary.csv")
        rank = pd.read_csv(ANALYSIS / f"{rec}_final_ranking.csv")
        dual_set = set()
        if rec == "NK2R":
            dual_set = set(rank[rank["both"] == True]["peptide_id"].tolist())
        for _, r in summ.iterrows():
            pid = r["peptide_id"]
            seq = clean(r["sequence"])
            known = pid in ("pep_known", "Substance_P", "Neurokinin_B")
            code = (f"{rec}-K" if known else
                    f"{rec}-{int(pid.replace('pep','')):02d}" if pid.startswith("pep") else pid)
            rows.append({
                "seq": seq,
                "group": f"ours_{rec}",
                "label": code,
                "rank": np.nan,
                "length": len(seq),
                "ec50": np.nan,
                "name": pid,
                "is_top10": False,
                "is_known": known,
                "is_dual": pid in dual_set,
                "is_own_team": False,
                "in_leaderboard": seq in lb_seqs,   # 该候选同时出现在排行榜前1000
                "team": "",
            })
    out = pd.DataFrame(rows)
    n_ours = int((~out["group"].eq("leaderboard")).sum())
    print(f"[OK] ours: {n_ours} 条 (其中与排行榜重合 {int(out['in_leaderboard'].fillna(False).sum())} 条)")

    # ---- 3. 去重: 我方候选优先保留 (便于高亮), 竞赛条目标记 in_leaderboard ----
    before = len(out)
    # 使我方候选排在前面 -> 去重时保留我方版本
    out["_pri"] = out["group"].eq("leaderboard").astype(int)
    out = out.sort_values("_pri", kind="stable").drop_duplicates("seq", keep="first")
    out = out.drop(columns="_pri").reset_index(drop=True)
    print(f"[OK] 去重: {before} -> {len(out)} (删除 {before-len(out)})")

    # ---- 4. 统计 ----
    print(f"\n序列长度分布: {out['length'].describe()[['min','25%','50%','75%','max']].to_dict()}")
    print(f"含 X 的序列: {(out['seq'].str.contains('X')).sum()}")
    print(f"\n分组统计:\n{out['group'].value_counts().to_string()}")

    out.to_csv(ANALYSIS / "umap_input.csv", index=False)
    print(f"\n[OK] 保存: {ANALYSIS}/umap_input.csv")


if __name__ == "__main__":
    main()
