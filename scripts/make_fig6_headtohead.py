#!/usr/bin/env python3
"""
Figure 6 (Main): 竞赛 Top10 vs 本研究 — 头对头对比总览

设计要点:
- a: 达标率对比 (横向分组条形图, 避免 x 轴标签重叠)
- b: 计分卡热图 (竞赛 Top10 + 我方最优 8 条, 指标归一化着色)
- c: z 值分布对比 (箱线+散点, 替代拥挤的 43 行排名)
- d: 长度分组表现 (图例移至图外, 避免遮挡)

配色: 竞赛=橙 #E28E2C, 本研究=蓝 #0F4D92, 语义色 绿=优/红=劣
"""
import matplotlib as mpl
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
import numpy as np
from pathlib import Path

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
# ---------------------------------------------------------------------------

# ---- palette ----
ORANGE = "#E28E2C"      # competition
ORANGE_L = "#F2C48A"
BLUE = PALETTE["blue_main"]        # this study
BLUE_L = "#7FA8D0"
GREEN = "#4C9F70"       # good / pass
RED = "#C0504D"         # poor / fail
GRAY = "#8A8A8A"
LGRAY = "#DCDCDC"
TEAL = PALETTE["teal"]

KNOWN_POS = 771.13

# heatmap: 红(差) - 白(中) - 绿(优) diverging
CMAP_PERF = LinearSegmentedColormap.from_list("perf", [RED, "#F5F0E8", GREEN])


def load_data():
    comp = pd.read_csv(ANALYSIS / "COMP_TOP10_summary.csv")
    ours = pd.read_csv(ANALYSIS / "NK2R_final_ranking.csv")
    ctl = pd.read_csv(ANALYSIS / "length_control.csv")

    ours["code"] = ["NK2R-K" if p == "pep_known" else f"NK2R-{int(p.replace('pep','')):02d}"
                    for p in ours["peptide_id"]]
    ours["length"] = 10
    stats = ctl.groupby("len")["pos"].agg(["mean", "std"])
    ours["z_vs_length"] = (ours["pos_score"] - stats.loc[10, "mean"]) / stats.loc[10, "std"]
    return comp, ours, ctl


def panel_a(ax, comp, ours, oo):
    """达标率对比 — 横向条形图"""
    metrics = [
        ("Dual-optimal\n(ΔPos<0 & frac_safer>0.5)", int(comp["both"].sum()) / 10 * 100,
         int(oo["both"].sum()) / len(oo) * 100),
        ("Off-target safe\n(frac_safer>0.5)", (comp["frac_safer"] > 0.5).sum() / 10 * 100,
         (oo["frac_safer"] > 0.5).sum() / len(oo) * 100),
        ("Stronger than known\n(ΔPos<0)", (comp["ΔPos"] < 0).sum() / 10 * 100,
         (ours["ΔPos"] < 0).sum() / 32 * 100),
        ("Beats random control\n(z<0)", (comp["z_vs_length"] < 0).sum() / 10 * 100,
         (ours["z_vs_length"] < 0).sum() / 32 * 100),
    ]
    y = np.arange(len(metrics))
    h = 0.34
    v1 = [m[1] for m in metrics]
    v2 = [m[2] for m in metrics]

    ax.barh(y + h/2, v1, h, color=ORANGE, edgecolor="white", linewidth=0.5,
            label="Competition Top10 (n=10)")
    ax.barh(y - h/2, v2, h, color=BLUE, edgecolor="white", linewidth=0.5,
            label="This study")

    for i in range(len(metrics)):
        ax.annotate(f"{v1[i]:.0f}", (v1[i] + 1.5, i + h/2), va="center", fontsize=5.5,
                    color=ORANGE, fontweight="bold")
        ax.annotate(f"{v2[i]:.0f}", (v2[i] + 1.5, i - h/2), va="center", fontsize=5.5,
                    color=BLUE, fontweight="bold")

    ax.set_yticks(y)
    ax.set_yticklabels([m[0] for m in metrics], fontsize=5.5)
    ax.set_xlim(0, 118)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Pass rate (%)")
    ax.axvline(50, color=LGRAY, linewidth=0.6, linestyle=":", zorder=0)
    ax.legend(loc="lower right", fontsize=5.5, handlelength=1.2, borderpad=0.3)
    ax.set_title("a  NK2R: criterion pass rates", loc="left", fontweight="bold", fontsize=8)
    ax.tick_params(length=2)


def panel_b(ax, comp, ours):
    """计分卡热图 — 竞赛 Top10 + 我方双优候选 (Top8, 数据完整)"""
    comp_s = comp.sort_values("z_vs_length")
    # 我方: 双优候选按 ΔPos 排序取前 8 (均为完整数据, 代表实际推荐)
    ours_pool = ours[(ours["code"] != "NK2R-K") & (ours["both"] == True)].copy()
    ours_pool["dpos"] = ours_pool["pos_score"] - KNOWN_POS
    ours_s = ours_pool.sort_values("dpos").head(8)

    n_comp, n_ours = len(comp_s), len(ours_s)
    codes = comp_s["code"].tolist() + ours_s["code"].tolist()

    # metric columns: (label, values_comp, values_ours, higher_is_better, is_binary)
    ec50 = comp_s["EC50_nM"].values.astype(float)
    ec50_log = np.log10(ec50)
    ec50_norm = 1 - (ec50_log - np.nanmin(ec50_log)) / (np.nanmax(ec50_log) - np.nanmin(ec50_log))

    cols = []
    # EC50: only competition has data
    cols.append(("EC50", np.concatenate([ec50_norm, np.full(n_ours, np.nan)]),
                 np.concatenate([ec50, np.full(n_ours, np.nan)]), "num"))
    # z-score
    zc, zo = comp_s["z_vs_length"].values, ours_s["z_vs_length"].values
    z_all = np.concatenate([zc, zo])
    z_norm = 1 - (z_all - z_all.min()) / (z_all.max() - z_all.min())
    cols.append(("z-score", z_norm, z_all, "num"))
    # ΔPos
    dc, do = comp_s["ΔPos"].values, (ours_s["pos_score"] - KNOWN_POS).values
    d_all = np.concatenate([dc, do])
    d_norm = 1 - (d_all - d_all.min()) / (d_all.max() - d_all.min())
    cols.append(("ΔPos", d_norm, d_all, "num"))
    # frac_safer
    fc = comp_s["frac_safer"].values
    fo = ours_s["frac_safer"].values
    f_all = np.concatenate([fc, fo])
    f_norm = (f_all - np.nanmin(f_all)) / (np.nanmax(f_all) - np.nanmin(f_all))
    cols.append(("frac_safer", f_norm, f_all, "num"))
    # dual
    bc = comp_s["both"].values.astype(float)
    bo = ours_s["both"].values.astype(float)
    b_all = np.concatenate([bc, bo])
    cols.append(("Dual", b_all, b_all, "bin"))

    mat = np.vstack([c[1] for c in cols]).T  # n_pep × n_metric

    # masked array: NaN -> gray
    masked = np.ma.masked_invalid(mat)
    cmap = CMAP_PERF.copy()
    cmap.set_bad(color="#EFEFEF")

    im = ax.imshow(masked, cmap=cmap, aspect="auto", vmin=0, vmax=1)

    # annotations — text color chosen for contrast against the cell colour
    def contrast_text(norm_val):
        if np.isnan(norm_val):
            return DGRAY
        # dark cells (norm near 0 red / near 1 green) need light text
        return "white" if (norm_val < 0.22 or norm_val > 0.78) else DGRAY

    for i in range(len(codes)):
        for j, (lab, norm_vals, raw, kind) in enumerate(cols):
            v = raw[i]
            nv = norm_vals[i]
            if np.isnan(v):
                txt = "n/a"
            elif kind == "bin":
                txt = "Y" if v == 1 else "N"
            elif lab == "EC50":
                txt = f"{v:.0f}"
            elif lab == "frac_safer":
                txt = f"{v:.2f}"
            else:
                txt = f"{v:+.1f}" if abs(v) < 10 else f"{v:+.0f}"
            ax.annotate(txt, (j, i), ha="center", va="center", fontsize=4.5,
                        color=contrast_text(nv), fontweight="bold" if kind == "bin" else "normal")

    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([c[0] for c in cols], fontsize=5.5, rotation=0)
    ax.set_yticks(range(len(codes)))
    ax.set_yticklabels(codes, fontsize=5)
    # color-code y labels by group
    for i, t in enumerate(ax.get_yticklabels()):
        t.set_color(ORANGE if i < n_comp else BLUE)
    ax.tick_params(length=0)

    # group divider + labels
    ax.axhline(n_comp - 0.5, color="black", linewidth=1.0)
    ax.annotate("Competition\nTop10", (-1.35, (n_comp - 1) / 2), ha="center", va="center",
                fontsize=5.5, color=ORANGE, fontweight="bold", rotation=90)
    ax.annotate("This study\n(dual-optimal)", (-1.35, n_comp + (n_ours - 1) / 2), ha="center",
                va="center", fontsize=5.5, color=BLUE, fontweight="bold", rotation=90)
    for s in ax.spines.values():
        s.set_visible(False)

    # colorbar
    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.028, pad=0.02, ticks=[0, 1])
    cbar.set_ticklabels(["poor", "good"])
    cbar.ax.tick_params(labelsize=5, length=0)
    cbar.outline.set_visible(False)

    ax.set_title("b  NK2R: scorecard (sorted by z-score)", loc="left", fontweight="bold", fontsize=8)


def panel_c(ax, comp, ours):
    """z 值分布对比 — 箱线 + 散点"""
    zc = comp["z_vs_length"].dropna().values
    zo = ours[ours["code"] != "NK2R-K"]["z_vs_length"].dropna().values

    bp = ax.boxplot([zc, zo], positions=[0, 1], widths=0.42, patch_artist=True,
                    showfliers=False, medianprops=dict(color=DGRAY, linewidth=1.0),
                    whiskerprops=dict(color=GRAY, linewidth=0.8),
                    capprops=dict(color=GRAY, linewidth=0.8))
    for box, fc, ec in zip(bp["boxes"], [ORANGE_L, BLUE_L], [ORANGE, BLUE]):
        box.set_facecolor(fc)
        box.set_alpha(0.45)
        box.set_edgecolor(ec)
        box.set_linewidth(1.0)

    rng = np.random.default_rng(0)
    jit = lambda n: rng.normal(0, 0.055, n)
    ax.scatter(np.zeros(len(zc)) + jit(len(zc)), zc, s=22, color=ORANGE,
               edgecolor="white", linewidth=0.5, zorder=3)
    ax.scatter(np.ones(len(zo)) + jit(len(zo)), zo, s=18, color=BLUE,
               edgecolor="white", linewidth=0.4, alpha=0.85, zorder=3)

    # 零线只画到数据区右侧, 给说明文字留位置 (避免文字压线)
    ax.plot([-0.5, 1.42], [0, 0], color=GRAY, linestyle="--", linewidth=0.8, zorder=1)
    ax.annotate("random-\ncontrol\nlevel", xy=(1.45, 0), xytext=(1.47, 0),
                fontsize=4.5, color=GRAY, va="center", ha="left")

    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        [f"Competition Top10\nn={len(zc)} · median {np.median(zc):+.2f}",
         f"This study\nn={len(zo)} · median {np.median(zo):+.2f}"], fontsize=5.5)
    ax.set_xlim(-0.5, 1.72)
    ax.set_ylabel("z-score vs length-matched random controls")
    # 中位数写进轴标签区, 避免与散点/图例打架
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax + 1.1)
    # 中位数并入 x 轴标签, 把顶部空间让给显著性括号 (避免两者叠压)
    ax.set_title("c  NK2R: binding evidence strength", loc="left", fontweight="bold", fontsize=8)
    # 显著性: 竞赛 Top10 vs 本研究 32 条候选 (双侧 Mann-Whitney U)
    _lo, _hi = ax.get_ylim()
    ax.set_ylim(_lo, _hi + (_hi - _lo) * 0.16)
    _p = mwu_p(zc, zo)
    sig_bracket(ax, 0, 1, _hi + (_hi - _lo) * 0.07, _p,
                label=f"n={len(zc)} vs n={len(zo)}")
    ax.tick_params(length=2)


def panel_d(ax, comp, oo, ours):
    """长度分组表现 — 图例移到图外"""
    c6 = comp[comp["length"] == 6]
    c9 = comp[comp["length"] == 9]
    c15 = comp[comp["length"] == 15]

    # n 已由条形上方的 "x/n" 标注体现, 轴标签保持两行避免与图例碰撞
    groups = [
        ("Competition\n6-mer", c6, ORANGE, len(c6)),
        ("Competition\n9-mer", c9, ORANGE, len(c9)),
        ("Competition\n15-mer", c15, ORANGE, len(c15)),
        ("This study\n10-mer", oo, BLUE, len(oo)),
    ]
    x = np.arange(len(groups))
    w = 0.26

    p_pos, p_safe, p_both = [], [], []
    for _, g, _, n in groups:
        p_pos.append((g["ΔPos"] < 0).sum() / n * 100)
        p_safe.append((g["frac_safer"] > 0.5).sum() / n * 100)
        p_both.append((g["both"] == True).sum() / n * 100)

    ax.bar(x - w, p_pos, w, color="#A8C8A0", edgecolor="white", linewidth=0.5,
           label="Delta-on < 0")
    ax.bar(x, p_safe, w, color=TEAL, edgecolor="white", linewidth=0.5,
           label="Fraction > 0.50")
    ax.bar(x + w, p_both, w, color=GREEN, edgecolor="white", linewidth=0.5,
           label="Primary rule")

    for i, (_, g, _, n) in enumerate(groups):
        ax.annotate(f"{(g['ΔPos']<0).sum()}/{n}", (i - w, p_pos[i] + 3), ha="center", fontsize=4.5)
        ax.annotate(f"{(g['frac_safer']>0.5).sum()}/{n}", (i, p_safe[i] + 3), ha="center", fontsize=4.5)
        ax.annotate(f"{(g['both']==True).sum()}/{n}", (i + w, p_both[i] + 3), ha="center", fontsize=4.5)

    ax.set_xticks(x)
    ax.set_xticklabels([g[0] for g in groups], fontsize=5.5)
    ax.set_ylim(0, 118)          # headroom for value annotations only
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("Pass rate (%)")
    # legend below the axes (no overlap with bars)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3, fontsize=5,
              handlelength=1.0, columnspacing=1.2, handletextpad=0.4)
    ax.set_title("d  NK2R: performance by length group", loc="left", fontweight="bold", fontsize=8)
    ax.tick_params(length=2)


def fig6():
    comp, ours, ctl = load_data()
    oo = ours.dropna(subset=["frac_safer"])

    fig = plt.figure(figsize=(183/25.4, 150/25.4))
    gs = fig.add_gridspec(2, 2, wspace=0.42, hspace=0.55,
                          height_ratios=[1.25, 1], width_ratios=[1, 1.05])

    panel_a(fig.add_subplot(gs[0, 0]), comp, ours, oo)
    panel_b(fig.add_subplot(gs[0, 1]), comp, ours)
    panel_c(fig.add_subplot(gs[1, 0]), comp, ours)
    panel_d(fig.add_subplot(gs[1, 1]), comp, oo, ours)

    fig.savefig(FIGDIR / "Fig6_head_to_head.svg", bbox_inches="tight")
    fig.savefig(FIGDIR / "Fig6_head_to_head.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / "Fig6_head_to_head.tiff", dpi=600, bbox_inches="tight")
    print("[OK] Fig6 saved (4 panels, redesigned)")


if __name__ == "__main__":
    fig6()
