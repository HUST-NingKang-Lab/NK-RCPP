#!/usr/bin/env python3
"""
Figure 3 (Main): 候选肽的脱靶选择性 — 4 panels
- a: 已知肽 vs 最优候选在全部脱靶蛋白上的配对分数散点 (3 receptors)
- b: 全部候选的 frac_safer 排名 (barh, >50% 脱靶蛋白结合更弱)
- c: 每受体 ΔNeg 中位数分布 (boxplot)
- d: 结构图: Substance P vs NK1R-21 在脱靶蛋白 3rze (H1 histamine receptor) 上的对接姿态
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
RAW = Path(os.environ.get("NK_RCPP_RAW", BASE / "raw"))
DOCKING = RAW / "peptide_docking"
RAPIDOCK = DOCKING / "RAPiDock"
RESULTS_DIR = RAPIDOCK / "results"
# ---------------------------------------------------------------------------
apply_style()

BLUE = PALETTE["blue_main"]
TEAL = PALETTE["teal"]
VIOLET = PALETTE["violet"]
GREEN = "#8BCF8B"
DGRAY = "#4D4D4D"    # 深灰, 用于强调文字
LIGHT = "#CFCECE"
RED = "#B64342"


RECEPTORS = {
    "NK1R": {"neg_dir": "NK1R_neg", "known_seq": "RPKPQQFFGLM",
             "best": "SPRRQQCFRL", "best_code": "NK1R-21", "color": BLUE},
    "NK2R": {"neg_dir": "NK2R_neg", "known_seq": "HKTDSFVGLM",
             "best": "KKNRQSCGGM", "best_code": "NK2R-10", "color": TEAL},
    "NK3R": {"neg_dir": "NK3R_neg", "known_seq": "DMHDFFVGLM",
             "best": "DEPTQQCIGF", "best_code": "NK3R-21", "color": VIOLET},
}
KNOWN_NAMES = {"NK1R": "Substance P", "NK2R": "NKA-derived", "NK3R": "Neurokinin B"}
BEST_PEP_ID = {"NK1R": "pep21", "NK2R": "pep10", "NK3R": "pep21"}

def read_best_score(score_csv):
    try:
        return float(pd.read_csv(score_csv)["ref2015score"].min())
    except Exception:
        return np.nan

def get_seq_scores(receptor):
    neg_dir = RESULTS_DIR / RECEPTORS[receptor]["neg_dir"]
    label_scores = {}
    for score_file in neg_dir.glob("*/**/ref2015_score.csv"):
        complex_name = score_file.parent.name
        if "_peptide" in complex_name:
            label = complex_name.split("_peptide")[-1]
        elif "_known" in complex_name:
            label = "known"
        else:
            continue
        protein = complex_name.rsplit("_", 1)[0]
        label_scores.setdefault(label, {})[protein] = read_best_score(score_file)
    neg_inputs = RESULTS_DIR.parent / f"{receptor.lower()}_neg_inputs"
    by_seq = {}
    for label, scores in label_scores.items():
        f = neg_inputs / ("known_peptide_docking_input.csv" if label == "known"
                          else f"peptide{label}_docking_input.csv")
        if f.exists():
            seq = str(pd.read_csv(f)["peptide_description"].iloc[0]).strip()
            by_seq.setdefault(seq, {}).update(scores)
    return by_seq

def fig3():
    fig = plt.figure(figsize=(183/25.4, 130/25.4))
    gs = fig.add_gridspec(2, 3, wspace=0.42, hspace=0.75,
                          width_ratios=[1, 1, 1.35])

    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])
    axc = fig.add_subplot(gs[0, 2])

    # --- a: mini paired scatter per receptor (inset style) ---
    for i, (rec, info) in enumerate(RECEPTORS.items()):
        seq_scores = get_seq_scores(rec)
        known = seq_scores[info["known_seq"]]
        best = seq_scores[info["best"]]
        common = sorted(set(known) & set(best))
        x = np.array([known[p] for p in common])
        y = np.array([best[p] for p in common])
        lim = [min(x.min(), y.min()) - 100, max(x.max(), y.max()) + 100]

        if i == 0:
            axa.scatter(x, y, s=8, color=info["color"], alpha=0.7, edgecolor="none")
            axa.plot(lim, lim, color=GRAY, linestyle="--", linewidth=0.6)
            axa.set_xlim(lim); axa.set_ylim(lim)
            axa.set_title(f"a — {rec} ({info['best_code']} vs control)",
                          loc="left", fontweight="bold", fontsize=7)
            axa.set_xlabel("Receptor-matched control score"); axa.set_ylabel("Candidate score")
        elif i == 1:
            axb.scatter(x, y, s=8, color=info["color"], alpha=0.7, edgecolor="none")
            axb.plot(lim, lim, color=GRAY, linestyle="--", linewidth=0.6)
            axb.set_xlim(lim); axb.set_ylim(lim)
            axb.set_title(f"b — {rec} ({info['best_code']} vs control)",
                          loc="left", fontweight="bold", fontsize=7)
            axb.set_xlabel("Receptor-matched control score"); axb.set_ylabel("Candidate score")
        else:
            axc.scatter(x, y, s=8, color=info["color"], alpha=0.7, edgecolor="none")
            axc.plot(lim, lim, color=GRAY, linestyle="--", linewidth=0.6)
            axc.set_xlim(lim); axc.set_ylim(lim)
            axc.set_title(f"c — {rec} ({info['best_code']} vs control)",
                          loc="left", fontweight="bold", fontsize=7)
            axc.set_xlabel("Receptor-matched control score"); axc.set_ylabel("Candidate score")
        for p_ax in (axa, axb, axc):
            p_ax.tick_params(labelsize=5.5)
        deltas = y - x
        n_safer = (deltas > 0).sum()
        # 数据集中在左下角(候选分数低), 故统计文字放右下空白角 + 白底
        (axa if i == 0 else axb if i == 1 else axc).text(
            0.97, 0.05, f"n={len(common)}\nΔmed {np.median(deltas):+.0f}\nsafer {n_safer}/{len(deltas)}",
            transform=(axa if i == 0 else axb if i == 1 else axc).transAxes,
            fontsize=5, va="bottom", ha="right",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.85))

    # ---------- d (hero): frac_safer ranking across all candidates ----------
    axd2 = fig.add_subplot(gs[1, 0])
    axd2_2 = fig.add_subplot(gs[1, 1])
    ax4 = fig.add_subplot(gs[1, 2])
    # 全部候选共 82 条 (有配对数据者), 逐条标注必然挤成一团
    # -> 只展示通过全文筛选标准的双优候选: both = (ΔPos<0 且 frac_safer>0.5)
    all_rows = []
    for rec in RECEPTORS:
        rk = pd.read_csv(ANALYSIS / f"{rec}_final_ranking.csv").dropna(subset=["frac_safer"])
        for _, r in rk.iterrows():
            pid = str(r["peptide_id"])
            code = (f"{rec}-K" if pid in ["Substance_P", "pep_known", "Neurokinin_B"]
                    else f"{rec}-{int(pid.replace('pep','')):02d}")
            all_rows.append({"rec": rec, "code": code, "frac": r["frac_safer"],
                             "delta_median": r["delta_median"], "both": r["both"]})
    ard = pd.DataFrame(all_rows)
    n_all = len(ard)
    ard = ard[ard["both"] == True].sort_values(["rec", "frac"], ascending=[True, False])
    colors = [RECEPTORS[r]["color"] for r in ard["rec"]]
    y = np.arange(len(ard))
    axd2.barh(y, ard["frac"], color=colors, height=0.72, edgecolor="none")
    axd2.plot([0.5, 0.5], [-0.6, len(ard) - 0.5], color=GRAY, linestyle="--", linewidth=0.8)
    axd2.text(0.485, -1.05, "0.5", fontsize=4.5, color=GRAY, ha="right", va="center")
    axd2.set_yticks(y)
    axd2.set_yticklabels(ard["code"], fontsize=5.5)
    axd2.set_ylim(len(ard) - 0.5, -1.5)      # 顶部留白给阈值标注
    axd2.set_xlim(0, 1.0)
    axd2.set_xlabel("Fraction with higher structural-neighbour score")
    # 标题保持短, 避免撞到右侧 e 面板标题; 计数放到图内右上空白处
    axd2.set_title("d — primary-rule candidates", loc="left", fontweight="bold", fontsize=8)
    axd2.text(0.985, -1.05, f"{len(ard)} of {n_all} candidates",
              transform=axd2.get_yaxis_transform(), ha="right", va="center",
              fontsize=5, color=DGRAY)
    axd2.tick_params(labelsize=5.5, length=2)

    # ---------- d2: ΔNeg median distribution by receptor ----------
    med_rows = []
    for rec in RECEPTORS:
        df = pd.read_csv(ANALYSIS / f"{rec}_pairwise_selectivity.csv").dropna(subset=["delta_median"])
        med_rows.append(df["delta_median"].tolist())
    bp = axd2_2.boxplot(med_rows, showfliers=False, widths=0.4, patch_artist=True,
                        medianprops=dict(color="black", linewidth=0.8))
    for i, (box, rec) in enumerate(zip(bp["boxes"], RECEPTORS)):
        box.set_facecolor("white")
        box.set_edgecolor(RECEPTORS[rec]["color"])
        box.set_linewidth(0.8)
    axd2_2.axhline(0, color=GRAY, linestyle="--", linewidth=0.7)
    axd2_2.set_xticks(range(1, 4))
    axd2_2.set_xticklabels(list(RECEPTORS.keys()), fontsize=6)
    axd2_2.set_ylabel("Mean matched candidate − control REF15 score")
    axd2_2.set_title("e — ΔNeg distribution", loc="left", fontweight="bold", fontsize=8)
    axd2_2.tick_params(labelsize=5.5, length=2)

    # ---------- f: binding mode on off-target 3rze ----------
    draw_offtarget_mode(ax4)
    ax4.set_title("f — poses on structural neighbour H1R (3rze)",
                  loc="left", fontweight="bold", fontsize=7)

    fig.savefig(FIGDIR / "Fig3_selectivity.svg", bbox_inches="tight")
    fig.savefig(FIGDIR / "Fig3_selectivity.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / "Fig3_selectivity.tiff", dpi=600, bbox_inches="tight")
    print("[OK] Fig3 saved (4 panels)")

def draw_offtarget_mode(ax):
    """Compare Substance P vs NK1R-21 (pep21) poses on off-target H1R (3rze)."""
    neg_dir = RESULTS_DIR / "NK1R_neg"
    known_pdb = neg_dir / "known_peptide_docking_input" / "3rze_A_known" / "rank1_ref2015.pdb"
    best_pdb = neg_dir / "peptide21_docking_input" / "3rze_A_peptide21" / "rank1_ref2015.pdb"

    def ca_trace(pdb_file):
        coords, resnames = [], []
        for line in open(pdb_file):
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                resnames.append(line[17:20].strip())
        return np.array(coords), resnames

    if not known_pdb.exists() or not best_pdb.exists():
        ax.text(0.5, 0.5, "poses unavailable", ha="center", va="center", fontsize=6)
        ax.set_xticks([]); ax.set_yticks([])
        return

    kc, kr = ca_trace(known_pdb)
    bc, br = ca_trace(best_pdb)
    # project to YZ
    ax.plot(kc[:, 1], kc[:, 2], "-o", color=GRAY, linewidth=1.2, markersize=3,
            markerfacecolor="white", markeredgewidth=0.6, label="Substance P")
    ax.plot(bc[:, 1], bc[:, 2], "-o", color=RED, linewidth=1.2, markersize=3,
            markerfacecolor="white", markeredgewidth=0.6, label="NK1R-21")
    # 端基标注加白底, 避免压在轨迹线上
    _bb = dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.8)
    for _c, _col in [(kc, GRAY), (bc, RED)]:
        ax.annotate("N", (_c[0, 1], _c[0, 2]), fontsize=4.5, color=_col, xytext=(0, -8),
                    textcoords="offset points", ha="center", zorder=8, bbox=_bb)
        ax.annotate("C", (_c[-1, 1], _c[-1, 2]), fontsize=4.5, color=_col, xytext=(0, -8),
                    textcoords="offset points", ha="center", zorder=8, bbox=_bb)
    ax.set_xticks([]); ax.set_yticks([])
    # 底部留白, 避免端基标注(尤其 N 端)越出坐标轴
    _ylo, _yhi = ax.get_ylim()
    ax.set_ylim(_ylo - (_yhi - _ylo) * 0.08, _yhi)
    # 图例移到坐标轴外下方, 不再压住轨迹
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=2,
              fontsize=5.5, columnspacing=1.5, handletextpad=0.5)

if __name__ == "__main__":
    fig3()
