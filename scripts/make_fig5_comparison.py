#!/usr/bin/env python3
"""
Figure 5 (Main): 竞赛 Top10 肽 vs 本研究的 NK2R 候选 — 对接性能与选择性对比
- a: 长度对照校正曲线 (ref2015 分数 vs 肽长, 随机组成匹配基线 + 竞赛肽 + 本研究)
- b: 正向 ref2015 分数分布 (三组)
- c: 实验 EC50 vs 长度校正后 z 值 (检验对接打分的预测力)
- d: ΔPos vs frac_safer 象限 (结合增益 × 脱靶安全)
- e: frac_safer 分布对比
- f: 结合姿态对比 (竞赛冠军 vs 本研究最优 vs 已知肽, NK2R 口袋)
"""
import matplotlib as mpl
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats as st

import sys
sys.path.insert(0, str(Path(__file__).parent))
from fig_style import (PALETTE, RECEPTOR_COLORS, AMBER, GRAY, DGRAY,
                       apply_style, save_pub)
apply_style()
from fig_style import sig_bracket, mwu_p

# --- paths ---------------------------------------------------------------
# Repository-relative: the scripts run from anywhere inside a checkout.
# Raw RAPiDock run directories and the competition workbook are not shipped;
# set NK_RCPP_RAW / NK_RCPP_COMP_XLSX (or NK_RCPP_BASE) to point at them.
# See scripts/README.md.
BASE = Path(os.environ.get("NK_RCPP_BASE", Path(__file__).resolve().parent.parent))
ANALYSIS = BASE / "step5_analysis" / "output"
FIGDIR = BASE / "figures"
RAW = Path(os.environ.get("NK_RCPP_RAW", BASE / "raw"))
DOCKING = RAW / "peptide_docking"
RAPIDOCK = DOCKING / "RAPiDock"
RESULTS_DIR = RAPIDOCK / "results"
# ---------------------------------------------------------------------------

BLUE = PALETTE["blue_main"]
TEAL = PALETTE["teal"]
GREEN = "#8BCF8B"
RED = "#B64342"
ORANGE = "#E28E2C"
LIGHT = "#CFCECE"

RESULTS = RESULTS_DIR
KNOWN_POS = 771.13


def fig5():
    comp = pd.read_csv(ANALYSIS / "COMP_TOP10_summary.csv")
    ours = pd.read_csv(ANALYSIS / "NK2R_final_ranking.csv")
    ctl = pd.read_csv(ANALYSIS / "length_control.csv")
    ours["code"] = ["NK2R-K" if p == "pep_known" else f"NK2R-{int(p.replace('pep','')):02d}"
                    for p in ours["peptide_id"]]
    ours["length"] = 10
    ours["z_vs_length"] = ours["pos_score"].apply(
        lambda s: (s - ctl[ctl["len"] == 10]["pos"].mean()) / ctl[ctl["len"] == 10]["pos"].std())
    # final_ranking.csv 含已知肽(基线), 本研究的统计必须排除它 —— 它单独作为基准展示
    ours_cand = ours[ours["code"] != "NK2R-K"].copy()
    assert len(ours_cand) == 32, f"预期 32 条候选, 实得 {len(ours_cand)}"

    fig = plt.figure(figsize=(183/25.4, 175/25.4))
    gs = fig.add_gridspec(3, 2, wspace=0.42, hspace=0.62)

    # ---------- a: length calibration ----------
    axa = fig.add_subplot(gs[0, 0])
    lengths = sorted(ctl["len"].unique())
    for L in lengths:
        d = ctl[ctl["len"] == L]["pos"]
        axa.scatter(np.full(len(d), L) + np.random.default_rng(L).normal(0, 0.10, len(d)),
                    d, s=16, color=LIGHT, edgecolor=GRAY, linewidth=0.3, zorder=2)
    # mean ± sd band per length
    ms = ctl.groupby("len")["pos"].agg(["mean", "std"])
    axa.errorbar(ms.index, ms["mean"], yerr=ms["std"], fmt="o-", color=GRAY,
                 linewidth=1.2, markersize=5, capsize=3, zorder=3, label="Random controls\n(mean ± s.d.)")
    # competition peptides — 标签按 y 排序后强制最小间距, 消除互相重叠
    cxs = (comp["length"] + 0.30).values
    cys = comp["pos_score"].values
    clabs = [c.replace("COMP-", "C") for c in comp["code"]]
    ylo, yhi = 550, 2350
    gap = (yhi - ylo) * 0.052
    order = np.argsort(cys)
    placed = []
    for i in order:
        yy = cys[i]
        for py in placed:
            if abs(yy - py) < gap:
                yy = py + gap
        placed.append(yy)
        # 9-mer 组若标签向右会压到 x=10 的蓝色菱形, 故整组改为向左放置
        if comp["length"].values[i] == 9:
            tx, ha = cxs[i] - 0.30, "right"
        else:
            tx, ha = cxs[i] + 0.22, "left"
        axa.annotate(clabs[i], xy=(cxs[i], cys[i]), xytext=(tx, yy),
                     textcoords="data", fontsize=4.5, color=ORANGE,
                     va="center", ha=ha, zorder=6,
                     arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.45,
                                     shrinkA=0, shrinkB=2))
    axa.scatter(cxs, cys, s=26, color=ORANGE, marker="s",
                edgecolor="white", linewidth=0.4, zorder=5, label="Competition Top10")
    # ours: 32 条候选同为 10-mer, 逐点绘制会糊成一条竖带 -> 用均值±标准差单点表示
    o_mean, o_std = ours_cand["pos_score"].mean(), ours_cand["pos_score"].std()
    axa.errorbar([10.0], [o_mean], yerr=[o_std], fmt="D", color=BLUE, markersize=6,
                 capsize=3, linewidth=1.2, zorder=5,
                 label=f"This study (10-mers, n={len(ours_cand)})")
    axa.axhline(KNOWN_POS, color=RED, linestyle="--", linewidth=0.8)
    # 已知肽标注放到虚线左端空白处, 避开 x=15 的橙色点
    axa.text(5.35, KNOWN_POS - 95, "receptor-matched control", fontsize=5, color=RED,
             va="top", ha="left")
    axa.set_xticks([6, 9, 10, 15])
    axa.set_xlim(5.2, 16.6)
    axa.set_ylim(ylo, yhi)
    axa.set_xlabel("Peptide length (residues)")
    axa.set_ylabel("REF15 total model score (NK2R)")
    axa.set_title("a — NK2R: peptide length strongly biases raw docking score", loc="left",
                  fontweight="bold", fontsize=8)
    axa.legend(loc="upper left", fontsize=5.5)

    # ---------- b: positive score distributions ----------
    axb = fig.add_subplot(gs[0, 1])
    data = [comp["pos_score"].values, ours_cand["pos_score"].values, [KNOWN_POS]]
    bp = axb.boxplot(data, widths=0.5, patch_artist=True, showfliers=False,
                     medianprops=dict(color="black", linewidth=0.9))
    for box, fc, ec in zip(bp["boxes"], [ORANGE, BLUE, LIGHT], [ORANGE, BLUE, GRAY]):
        box.set_facecolor(fc); box.set_alpha(0.55); box.set_edgecolor(ec); box.set_linewidth(0.9)
    for i, d in enumerate(data, 1):
        jitter = np.random.default_rng(i).normal(0, 0.055, len(d))
        axb.scatter(np.full(len(d), i) + jitter, d, s=13, color=[ORANGE, BLUE, GRAY][i-1],
                    alpha=0.75, edgecolor="none", zorder=3)
    axb.axhline(KNOWN_POS, color=GRAY, linestyle="--", linewidth=0.8)
    axb.set_xticks([1, 2, 3])
    axb.set_xticklabels(["Competition\nTop10", "This study\n(32 cand.)", "Receptor-matched\ncontrol"], fontsize=6)
    axb.set_ylabel("REF15 total model score")
    axb.set_title("b — NK2R: raw positive docking scores", loc="left",
                  fontweight="bold", fontsize=8)
    # 显著性: 竞赛 Top10 vs 本研究候选 (双侧 Mann-Whitney U)
    _lo, _hi = axb.get_ylim()
    axb.set_ylim(_lo, _hi + (_hi - _lo) * 0.13)
    _p = mwu_p(comp["pos_score"], ours_cand["pos_score"])
    sig_bracket(axb, 1, 2, _hi + (_hi - _lo) * 0.055, _p,
                label=f"n={len(comp)} vs n={len(ours_cand)}")

    # ---------- c: EC50 vs length-corrected score ----------
    axc = fig.add_subplot(gs[1, 0])
    sub = comp.dropna(subset=["EC50_nM", "z_vs_length"])
    axc.scatter(sub["EC50_nM"], sub["z_vs_length"], s=30, color=ORANGE,
                edgecolor="white", linewidth=0.5, zorder=3)
    # 标签按 z 值排序后强制最小间距 (右下角 C05/C06/C08/C09 原本会叠在一起)
    _items = sorted(zip(sub["EC50_nM"].values, sub["z_vs_length"].values,
                        [c.replace("COMP-", "C") for c in sub["code"]]),
                    key=lambda t: t[1])
    _zlo, _zhi = -1.55, 1.75
    _gap = (_zhi - _zlo) * 0.085
    _placed = []
    for xv, yv, lb in _items:
        yy = yv
        for py in _placed:
            if abs(yy - py) < _gap:
                yy = py + _gap
        _placed.append(yy)
        axc.annotate(lb, xy=(xv, yv), xytext=(xv * 1.13, yy), textcoords="data",
                     fontsize=5, color=GRAY, va="center", ha="left", zorder=6,
                     arrowprops=dict(arrowstyle="-", color=GRAY, lw=0.45,
                                     shrinkA=0, shrinkB=2))
    axc.set_ylim(_zlo, _zhi)
    axc.axhline(0, color=GRAY, linestyle="--", linewidth=0.8)
    # NB: use axes-fraction coords — data coords on a log axis are unsafe before draw
    axc.annotate(" random-control level", xy=(0.02, 0.0), xycoords=("axes fraction", "data"),
                 fontsize=5, color=GRAY, va="bottom")
    axc.set_xscale("log")
    axc.set_xlabel("Experimental EC50 (nM, log scale)")
    axc.set_ylabel("z-score vs length-matched\nrandom controls")
    rho, pval = st.spearmanr(sub["EC50_nM"], sub["z_vs_length"])
    # 左下/右下均有数据点, 统计文字放左上空白区 (EC50 2–20 nM, z 0–1.5 区间无点) + 白底
    axc.annotate(f"Spearman ρ = {rho:+.2f}\n(P = {pval:.2f}, n={len(sub)})",
                 xy=(0.025, 0.97), xycoords="axes fraction", ha="left", va="top",
                 fontsize=5.5,
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.85))
    axc.set_title("c — NK2R: docking score vs experimental potency", loc="left", fontweight="bold", fontsize=8)

    # ---------- d: ΔPos vs frac_safer ----------
    axd = fig.add_subplot(gs[1, 1])
    axd.axhline(0.5, color=GRAY, linestyle="--", linewidth=0.7)
    axd.axvline(0, color=GRAY, linestyle="--", linewidth=0.7)
    axd.axvspan(-400, 0, color=GREEN, alpha=0.07)
    axd.axhspan(0.5, 1.0, color=GREEN, alpha=0.07)
    oc = comp.dropna(subset=["ΔPos", "frac_safer"])
    axd.scatter(oc["ΔPos"], oc["frac_safer"], s=26, color=ORANGE, marker="s",
                edgecolor="white", linewidth=0.5, label="Competition Top10", zorder=3)
    # 此面板橙色方块与蓝点交织, 逐点标号会压在蓝点上 -> 不标号(身份见 a、c 面板)
    oo = ours_cand.dropna(subset=["ΔPos", "frac_safer"])
    axd.scatter(oo["ΔPos"], oo["frac_safer"], s=20, color=BLUE, edgecolor="white",
                linewidth=0.4, alpha=0.85, label="This study", zorder=3)
    axd.set_xlabel("Delta-on (candidate − control)")
    axd.set_ylabel("Fraction with higher\nstructural-neighbour score")
    axd.legend(loc="lower right", fontsize=6)
    axd.set_title("d — NK2R: Delta-on vs structural-neighbour fraction", loc="left", fontweight="bold", fontsize=8)

    # ---------- e: frac_safer distribution ----------
    axe = fig.add_subplot(gs[2, 0])
    d1, d2 = oc["frac_safer"].values, oo["frac_safer"].values
    bp2 = axe.boxplot([d1, d2], widths=0.5, patch_artist=True, showfliers=False,
                      medianprops=dict(color="black", linewidth=0.9))
    for box, fc, ec in zip(bp2["boxes"], [ORANGE, BLUE], [ORANGE, BLUE]):
        box.set_facecolor(fc); box.set_alpha(0.55); box.set_edgecolor(ec); box.set_linewidth(0.9)
    for i, d in enumerate([d1, d2], 1):
        jitter = np.random.default_rng(i).normal(0, 0.06, len(d))
        axe.scatter(np.full(len(d), i) + jitter, d, s=14, color=[ORANGE, BLUE][i-1],
                    alpha=0.75, edgecolor="none", zorder=3)
    axe.axhline(0.5, color=GRAY, linestyle="--", linewidth=0.8)
    axe.set_xticks([1, 2])
    axe.set_xticklabels([f"Competition Top10\n(n={len(d1)})", f"This study\n(n={len(d2)})"], fontsize=6)
    axe.set_ylabel("Fraction with higher score")
    axe.set_title("e — NK2R: structural-neighbour fraction", loc="left",
                  fontweight="bold", fontsize=8)
    # 显著性: 竞赛 Top10 vs 本研究候选 (有聚合值者)
    _lo, _hi = axe.get_ylim()
    axe.set_ylim(_lo, _hi + (_hi - _lo) * 0.13)
    _p = mwu_p(oc["frac_safer"], oo["frac_safer"])
    sig_bracket(axe, 1, 2, _hi + (_hi - _lo) * 0.055, _p,
                label=f"n={len(oc)} vs n={len(oo)}")

    # ---------- f: binding pose overlay on NK2R ----------
    axf = fig.add_subplot(gs[2, 1])
    draw_pose_overlay(axf)
    axf.set_title("f — NK2R pocket: pose overlay (Cα traces)", loc="left", fontweight="bold", fontsize=8)

    fig.savefig(FIGDIR / "Fig5_competition_comparison.svg", bbox_inches="tight")
    fig.savefig(FIGDIR / "Fig5_competition_comparison.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / "Fig5_competition_comparison.tiff", dpi=600, bbox_inches="tight")
    print("[OK] Fig5 saved (6 panels)")


def draw_pose_overlay(ax):
    """Overlay Cα traces: competition champion (comp01), our best (pep27), known ligand."""
    def ca_trace(pdb_file):
        coords = []
        for line in open(pdb_file):
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        return np.array(coords)

    paths = {
        "Control (NKA-derived)": (RESULTS / "NK2R_peptide_len10/pep_known/rank1_ref2015.pdb", GRAY),
        "Best of this study (NK2R-27)": (RESULTS / "NK2R_peptide_len10/pep27/rank1_ref2015.pdb", BLUE),
        "Competition champion (COMP-01)": (RESULTS / "COMP_TOP10_NK2R/comp01/rank1_ref2015.pdb", ORANGE),
    }
    for label, (p, col) in paths.items():
        if not p.exists():
            continue
        c = ca_trace(p)
        if len(c) == 0:
            continue
        ax.plot(c[:, 1], c[:, 2], "-o", color=col, linewidth=1.2, markersize=3,
                markerfacecolor="white", markeredgewidth=0.5, label=label)
    ax.set_xticks([]); ax.set_yticks([])
    # 图例移到坐标轴外下方, 不再压住轨迹
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.03), ncol=2,
              fontsize=5.5, columnspacing=1.2, handletextpad=0.5)


if __name__ == "__main__":
    fig5()
