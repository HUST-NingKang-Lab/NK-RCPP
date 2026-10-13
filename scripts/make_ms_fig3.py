#!/usr/bin/env python3
"""
Manuscript Fig. 3 | NK-RCPP candidate-level structural-neighbour fraction rankings
and inherited matched score differences.

a–c  Evaluable candidates ranked within NK1R/NK2R/NK3R by the fraction of reported
     structural-neighbour comparisons in which the candidate REF15 total model score
     exceeded the receptor-matched control score.
     Receptor-coloured = primary joint rule; grey = not; star = strict >=0.75.
     Dashed 0.50, dotted amber 0.75.
     Titles report inherited comparison count / receptor-specific Table S3 row count.
d–f  Inherited mean matched candidate-minus-control REF15 score differences for the
     same candidates (same x order). Aggregates could not be reconstructed or linked
     to individual structures (pair-level rows absent); retrieved proteins were not
     validated pharmacological off-targets.

可读性改进 (相对初稿):
  - 图放宽至双栏 183 mm
  - x 轴候选代号隔位标注 (旋转 90°), 相邻标签间距翻倍, 不再互相叠压
  - 未标注的候选保留刻度线便于对位
"""
import matplotlib as mpl
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent))
from fig_style import (PALETTE, RECEPTOR_COLORS, AMBER, GRAY, DGRAY,
                       apply_style, save_pub)

# --- paths ---------------------------------------------------------------
# Repository-relative: the scripts run from anywhere inside a checkout.
# Raw RAPiDock run directories and the competition workbook are not shipped;
# set NK_RCPP_RAW / NK_RCPP_COMP_XLSX (or NK_RCPP_BASE) to point at them.
# See scripts/README.md.
BASE = Path(os.environ.get("NK_RCPP_BASE", Path(__file__).resolve().parent.parent))
ANALYSIS = BASE / "step5_analysis" / "output"
FIGDIR = BASE / "figures"
# ---------------------------------------------------------------------------
apply_style()

C = RECEPTOR_COLORS


CONTROL = {"NK1R": "Substance_P", "NK2R": "pep_known", "NK3R": "Neurokinin_B"}
TABLE_S3_ROWS = {"NK1R": 66, "NK2R": 27, "NK3R": 60}
RECS = ["NK1R", "NK2R", "NK3R"]


def load(rec):
    df = pd.read_csv(ANALYSIS / f"{rec}_final_ranking.csv")
    df = df[df["peptide_id"] != CONTROL[rec]].copy()
    df = df.dropna(subset=["frac_safer"])          # 仅 evaluable 候选
    df["code"] = [f"{rec}-{int(p.replace('pep','')):02d}" for p in df["peptide_id"]]
    df["primary"] = df["both"] == True
    df["strict"] = df["frac_safer"] >= 0.75
    df = df.sort_values("frac_safer", ascending=False).reset_index(drop=True)
    df["rank"] = np.arange(len(df))
    return df


def xaxis_labels(ax, codes, every=2):
    """隔位标注候选代号 (旋转 90°), 保持全部刻度线便于对位。"""
    n = len(codes)
    ax.set_xticks(np.arange(n))
    ax.set_xticklabels([codes[i] if i % every == 0 else "" for i in range(n)],
                       rotation=90, fontsize=5.5, va="top")
    ax.tick_params(axis="x", length=2, pad=2)


def panel_rank(ax, rec, df, letter):
    ax.axhline(0.50, color=DGRAY, linestyle="--", linewidth=0.8, zorder=1)
    ax.axhline(0.75, color=AMBER, linestyle=":", linewidth=1.0, zorder=1)

    ni = df[~df["primary"]]
    pr = df[df["primary"]]
    ax.scatter(ni["rank"], ni["frac_safer"], s=26, c=GRAY, zorder=3)
    ax.scatter(pr["rank"], pr["frac_safer"], s=30, c=C[rec], zorder=4)
    st = df[df["strict"]]
    if len(st):
        ax.scatter(st["rank"], st["frac_safer"], s=150, marker="*",
                   facecolors="none", edgecolors=DGRAY, linewidths=0.9, zorder=5)

    ax.set_ylim(0.28, 1.02)
    ax.set_yticks([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.set_xlim(-0.8, len(df) - 0.2)
    ax.set_ylabel("Fraction with higher\nstructural-neighbour score", fontsize=7)
    ax.tick_params(axis="y", labelsize=6, length=2)
    ax.set_title(f"{letter}  {rec}  (n={len(df)}; reported {TABLE_S3_ROWS[rec]}/{TABLE_S3_ROWS[rec]})",
                 loc="left", fontweight="bold", fontsize=8, color=C[rec])


def panel_delta(ax, rec, df, letter):
    ax.axhline(0, color=DGRAY, linestyle="--", linewidth=0.8, zorder=1)
    ni = df[~df["primary"]]
    pr = df[df["primary"]]
    ax.scatter(ni["rank"], ni["delta_mean"], s=26, c=GRAY, zorder=3)
    ax.scatter(pr["rank"], pr["delta_mean"], s=30, c=C[rec], zorder=4)
    ax.set_xlim(-0.8, len(df) - 0.2)
    ax.set_ylabel("Candidate − control\nREF15 total score", fontsize=7)
    ax.tick_params(axis="y", labelsize=6, length=2)
    ax.set_title(f"{letter}  Inherited matched mean delta", loc="left",
                 fontweight="bold", fontsize=8, color=DGRAY)


def main():
    data = [(r, load(r)) for r in RECS]

    fig = plt.figure(figsize=(183/25.4, 165/25.4))
    gs = fig.add_gridspec(2, 3, hspace=0.52, wspace=0.30)

    for i, (rec, df) in enumerate(data):
        ax = fig.add_subplot(gs[0, i])
        panel_rank(ax, rec, df, "abc"[i])
        xaxis_labels(ax, df["code"].tolist(), every=2)

        ax2 = fig.add_subplot(gs[1, i])
        panel_delta(ax2, rec, df, "def"[i])
        xaxis_labels(ax2, df["code"].tolist(), every=2)

    # 图例 (仅一次)
    from matplotlib.lines import Line2D
    leg = [
        Line2D([], [], marker="o", ls="", ms=4.5, mfc=C["NK1R"], mec="none",
               label="Primary rule"),
        Line2D([], [], marker="o", ls="", ms=4.5, mfc=GRAY, mec="none",
               label="Not retained"),
        Line2D([], [], marker="*", ls="", ms=9, mfc="none", mec=DGRAY, mew=0.9,
               label="Strict ≥0.75"),
        Line2D([], [], ls="--", color=DGRAY, lw=0.8, label="0.50"),
        Line2D([], [], ls=":", color=AMBER, lw=1.0, label="0.75"),
    ]
    fig.legend(handles=leg, loc="lower center", ncol=5, fontsize=6,
               bbox_to_anchor=(0.5, -0.015), columnspacing=1.4, handletextpad=0.5)

    fig.savefig(FIGDIR / "Fig3_audit.svg", bbox_inches="tight")
    fig.savefig(FIGDIR / "Fig3_audit.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / "Fig3_audit.tiff", dpi=600, bbox_inches="tight")
    print("[OK] Fig3_audit saved")
    for rec, df in data:
        print(f"   {rec}: n={len(df)} primary={int(df['primary'].sum())} strict={int(df['strict'].sum())}")


if __name__ == "__main__":
    main()
