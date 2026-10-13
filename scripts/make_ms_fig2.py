#!/usr/bin/env python3
"""
Manuscript Fig. 2 | Receptor-stratified audit of NK-RCPP on-target scores,
structural-neighbour aggregates and fixed prioritization rules.

a  Candidate-minus-control minimum REF15 total model scores (32 candidates per receptor)
b  Joint display of Delta-on vs fraction with higher structural-neighbour score,
   with count-proxy encoding, primary-rule outlines and aggregate-unavailable crosses
c  Counts: adequate count proxy / missing-low proxy / primary rule / strict >=0.75
d  Top primary-rule candidate per receptor (ordered by fraction then Delta-on)

术语严格遵循稿件:
  Delta-on             候选分数 − 受体匹配对照分数
  REF15 total model score   (非结合自由能)
  structural-neighbour  (非 validated pharmacological off-target)
  论文强调: fractions are inherited computational summaries, not biological safety probabilities.
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

# 与稿件 Fig 2 一致的受体配色

SEED = 20260801               # 稿件声明的固定随机种子
CONTROL = {"NK1R": "Substance_P", "NK2R": "pep_known", "NK3R": "Neurokinin_B"}
TABLE_S3_ROWS = {"NK1R": 66, "NK2R": 27, "NK3R": 60}
RECS = ["NK1R", "NK2R", "NK3R"]


def load(rec):
    """候选表 (排除受体匹配对照), 补 code / Delta-on / primary / strict / proxy."""
    df = pd.read_csv(ANALYSIS / f"{rec}_final_ranking.csv")
    df = df[df["peptide_id"] != CONTROL[rec]].copy()
    df["code"] = [f"{rec}-{int(p.replace('pep','')):02d}" for p in df["peptide_id"]]
    df["delta_on"] = df["pos_score"] - df.loc[df["peptide_id"] == CONTROL[rec], "pos_score"].iloc[0] \
        if (df["peptide_id"] == CONTROL[rec]).any() else df["ΔPos"]
    # ΔPos 已由上游算好(= 候选 − 对照), 直接采用
    df["delta_on"] = df["ΔPos"]
    df["primary"] = df["both"] == True
    df["strict"] = df["frac_safer"] >= 0.75
    df["proxy"] = df["n_common"] / TABLE_S3_ROWS[rec]
    return df


def panel_a(ax, rec, df):
    """Delta-on 分布 (strip plot)。彩色 = 满足 primary rule。"""
    rng = np.random.default_rng(SEED)
    d = df.dropna(subset=["delta_on"])
    ok = d[d["primary"]]
    no = d[~d["primary"]]
    ax.axvline(0, color=DGRAY, linestyle="--", linewidth=0.8, zorder=1)
    ax.scatter(no["delta_on"], rng.normal(0, 0.075, len(no)), s=22, c=GRAY,
               edgecolors="white", linewidths=0.4, zorder=2)
    ax.scatter(ok["delta_on"], rng.normal(0, 0.075, len(ok)), s=26, c=C[rec],
               edgecolors="white", linewidths=0.4, zorder=3)
    ax.set_ylim(-0.55, 0.55)
    ax.set_yticks([])
    ax.set_xlabel("Candidate − control\nREF15 total score")
    ax.tick_params(labelsize=6, length=2)
    lo, hi = d["delta_on"].min(), d["delta_on"].max()
    pad = (hi - lo) * 0.12 + 1
    ax.set_xlim(lo - pad, hi + pad)
    ax.annotate(f"n={len(d)}", xy=(0.02, 0.90), xycoords="axes fraction",
                fontsize=5.5, color=DGRAY, va="top")


def panel_b(ax, rec, df):
    """Delta-on × structural-neighbour fraction, 含 proxy 编码与 missing 叉号。"""
    d = df.dropna(subset=["frac_safer"])
    miss = df[df["frac_safer"].isna()]

    ax.axhline(0.50, color=DGRAY, linestyle="--", linewidth=0.8, zorder=1)
    ax.axhline(0.75, color=AMBER, linestyle=":", linewidth=1.0, zorder=1)
    ax.axvline(0, color=DGRAY, linestyle="--", linewidth=0.8, zorder=1)

    # 有聚合数据的候选: 实心 = proxy>=90% (本数据集全部 100%)
    good = d[d["proxy"] >= 0.90]
    bad = d[d["proxy"] < 0.90]
    ax.scatter(good["delta_on"], good["frac_safer"], s=26, c=C[rec],
               edgecolors="white", linewidths=0.4, zorder=3)
    if len(bad):
        ax.scatter(bad["delta_on"], bad["frac_safer"], s=26, facecolors="white",
                   edgecolors=C[rec], linewidths=0.8, zorder=3)
    # primary rule 候选加深色描边
    ok = d[d["primary"]]
    ax.scatter(ok["delta_on"], ok["frac_safer"], s=95, facecolors="none",
               edgecolors=DGRAY, linewidths=0.9, zorder=4)

    # 无聚合数据的候选: 底部叉号 (不给它们编造 fraction)
    if len(miss):
        yx = 0.075
        ax.scatter(miss["delta_on"], np.full(len(miss), yx), s=30, marker="x",
                   c=GRAY, linewidths=1.0, zorder=3)
        ax.annotate("× aggregate unavailable", xy=(0.02, 0.155), xycoords="axes fraction",
                    fontsize=5, color=DGRAY, va="bottom",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))
    else:
        # 说明文字加白底, 避免与 Δ=0 虚线交叠; house palette 的 neutral 作文字过浅, 改用 DGRAY
        ax.annotate("reported count proxy = 100%", xy=(0.02, 0.035),
                    xycoords="axes fraction", fontsize=5, color=DGRAY, va="bottom",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))

    ax.set_ylim(-0.01, 1.02)
    ax.set_yticks([0.0, 0.25, 0.50, 0.75, 1.00])
    ax.set_xlabel("Candidate − control total score")
    ax.tick_params(labelsize=6, length=2)
    lo, hi = df["delta_on"].min(), df["delta_on"].max()
    pad = (hi - lo) * 0.10 + 1
    ax.set_xlim(lo - pad, hi + pad)


def panel_c(ax, data):
    """计数条形: adequate proxy / missing-low / primary / strict。"""
    cats = ["Reported count\nproxy ≥90%", "Missing / low\ncount proxy",
            "Primary\nrule", "Strict\n≥0.75"]
    vals = {
        "Reported count\nproxy ≥90%": [int((d["proxy"] >= 0.90).sum()) for _, d in data],
        "Missing / low\ncount proxy": [int(d["proxy"].isna().sum()) for _, d in data],
        "Primary\nrule": [int(d["primary"].sum()) for _, d in data],
        "Strict\n≥0.75": [int(d["strict"].sum()) for _, d in data],
    }
    x = np.arange(len(cats))
    w = 0.26
    for i, rec in enumerate(RECS):
        off = (i - 1) * w
        ax.bar(x + off, [vals[c][i] for c in cats], w, color=C[rec],
               edgecolor="white", linewidth=0.5, label=rec)
        for j, c in enumerate(cats):
            v = vals[c][i]
            if v > 0:
                ax.annotate(str(v), (x[j] + off, v + 0.6), ha="center", fontsize=5.5,
                            color=DGRAY)
    ax.set_xticks(x)
    ax.set_xticklabels(cats, fontsize=6)
    ax.set_ylabel("Candidates")
    ax.set_ylim(0, 37)
    ax.tick_params(labelsize=6, length=2)
    ax.legend(loc="upper right", fontsize=6, ncol=3, columnspacing=0.9, handlelength=1.0)


def panel_d(ax, top):
    """各受体 primary-rule 最优候选 (按 fraction 再按 Delta-on 排序)。"""
    y = np.arange(len(top))
    ax.axvline(0.50, color=DGRAY, linestyle="--", linewidth=0.8, zorder=1)
    ax.axvline(0.75, color=AMBER, linestyle=":", linewidth=1.0, zorder=1)
    for i, (rec, row) in enumerate(top):
        ax.plot([0.50, row["frac_safer"]], [i, i], color=C[rec], linewidth=1.2, zorder=2)
        ax.scatter([row["frac_safer"]], [i], s=48, c=C[rec], edgecolors="white",
                   linewidths=0.6, zorder=3)
        ax.annotate(f"Δ={row['delta_on']:.1f}; {int(row['n_common'])}/{TABLE_S3_ROWS[rec]}",
                    (row["frac_safer"] + 0.012, i), va="center", fontsize=6, color=DGRAY)
    ax.set_yticks(y)
    ax.set_yticklabels([r["code"] for _, r in top], fontsize=7)
    ax.set_ylim(len(top) - 0.5, -0.5)
    ax.set_xlim(0.50, 1.02)
    ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.set_xlabel("Fraction with higher score")
    ax.tick_params(labelsize=6, length=2)


def main():
    data = [(r, load(r)) for r in RECS]

    fig = plt.figure(figsize=(183/25.4, 150/25.4))
    gs = fig.add_gridspec(3, 3, height_ratios=[0.85, 1.15, 0.85],
                          hspace=0.85, wspace=0.30)

    for i, (rec, df) in enumerate(data):
        axa = fig.add_subplot(gs[0, i])
        panel_a(axa, rec, df)
        axa.set_title(f"a  {rec}" if i == 0 else rec, loc="left" if i == 0 else "center",
                      fontweight="bold", fontsize=8,
                      color=DGRAY if i == 0 else C[rec])

        axb = fig.add_subplot(gs[1, i])
        panel_b(axb, rec, df)
        if i == 0:
            axb.set_ylabel("Fraction with higher\nstructural-neighbour score", fontsize=7)
        if i == 1:
            axb.annotate("b", xy=(0.02, 1.06), xycoords="axes fraction",
                         fontweight="bold", fontsize=9)

    axc = fig.add_subplot(gs[2, 0:2])
    panel_c(axc, data)
    axc.set_title("c", loc="left", fontweight="bold", fontsize=9)

    axd = fig.add_subplot(gs[2, 2])
    top = []
    for rec, df in data:
        cand = df[df["primary"]].dropna(subset=["frac_safer"])
        sel = cand.sort_values(["frac_safer", "delta_on"], ascending=[False, True]).iloc[0]
        top.append((rec, sel))
    panel_d(axd, top)
    axd.set_title("d  Top primary-rule candidate\n   (reported comparison counts)",
                  loc="left", fontweight="bold", fontsize=8, color=DGRAY)

    fig.savefig(FIGDIR / "Fig2_audit.svg", bbox_inches="tight")
    fig.savefig(FIGDIR / "Fig2_audit.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / "Fig2_audit.tiff", dpi=600, bbox_inches="tight")
    print("[OK] Fig2_audit saved")
    for rec, sel in top:
        print(f"   {sel['code']}: Δ={sel['delta_on']:.2f}, {int(sel['n_common'])}/{TABLE_S3_ROWS[rec]}, frac={sel['frac_safer']:.3f}")


if __name__ == "__main__":
    main()
