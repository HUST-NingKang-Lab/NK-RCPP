#!/usr/bin/env python3
"""
Figure 2 (Main): 三受体候选肽正向对接性能 — 4 panels
- a: 已知肽 vs 每个受体最优候选的 ref2015 分数 (bar)
- b: 全部候选 vs 已知肽的正向分数分布 (strip + known line)
- c: ΔPos vs ΔNeg 象限散点 (hero panel: 结合改善 vs 脱靶安全)
- d: 结合模式图: 代表性受体 (NK2R) 已知肽 vs 最优候选的对接姿态对比
"""
import matplotlib as mpl
import os
import matplotlib.pyplot as plt
import pandas as pd
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
RAW = Path(os.environ.get("NK_RCPP_RAW", BASE / "raw"))
DOCKING = RAW / "peptide_docking"
RAPIDOCK = DOCKING / "RAPiDock"
RESULTS_DIR = RAPIDOCK / "results"
# ---------------------------------------------------------------------------

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Liberation Sans', 'DejaVu Sans']
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 7
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['legend.frameon'] = False

BLUE = "#0F4D92"
BLUE_SOFT = "#3775BA"
TEAL = "#42949E"
VIOLET = "#9A4D8E"
GREEN = "#8BCF8B"
RED = "#B64342"
GOLD = "#E2A93B"
GRAY = "#767676"
LIGHT = "#CFCECE"
DGREEN = "#3F8F5F"   # 深绿, 用于浅绿底上的文字
DGRAY = "#4D4D4D"    # 深灰, 用于强调文字

FIGDIR.mkdir(exist_ok=True)

RECEPTOR_COLORS = {"NK1R": BLUE, "NK2R": TEAL, "NK3R": VIOLET}
KNOWN_POS = {"NK1R": 574.39, "NK2R": 771.13, "NK3R": 2240.54}
KNOWN_NAMES = {"NK1R": "Substance P", "NK2R": "NKA-derived", "NK3R": "Neurokinin B"}
BEST_PEP = {"NK1R": "pep18", "NK2R": "pep27", "NK3R": "pep21"}
BEST_CODE = {"NK1R": "NK1R-18", "NK2R": "NK2R-27", "NK3R": "NK3R-21"}
BEST_SEQ = {"NK1R": "SPRRQLGFGL", "NK2R": "PSNTQCCFGL", "NK3R": "DEPTQQCIGF"}

def load_final(receptor):
    df = pd.read_csv(ANALYSIS / f"{receptor}_final_ranking.csv")
    df = df[~df["peptide_id"].isin(["Substance_P", "pep_known", "Neurokinin_B"])]
    df = df.copy()
    df["ΔPos"] = pd.to_numeric(df["ΔPos"], errors="coerce")
    df["ΔNeg"] = pd.to_numeric(df["ΔNeg"], errors="coerce")
    # add code column
    df["code"] = [f"{receptor}-{int(pid.replace('pep','')):02d}" for pid in df["peptide_id"]]
    return df

def fig2():
    fig = plt.figure(figsize=(183/25.4, 130/25.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], wspace=0.4, hspace=0.65)
    receptors = ["NK1R", "NK2R", "NK3R"]

    # ---------- Panel a: known vs best candidate ----------
    ax = fig.add_subplot(gs[0, 0])
    bests = {}
    for rec in receptors:
        df = load_final(rec)
        bests[rec] = df.dropna(subset=["ΔPos"]).sort_values("ΔPos").iloc[0]["pos_score"]

    x = np.arange(len(receptors))
    width = 0.35
    known_vals = [KNOWN_POS[r] for r in receptors]
    best_vals = [bests[r] for r in receptors]

    ax.bar(x - width/2, known_vals, width, color=LIGHT, edgecolor="none", label="Known ligand")
    ax.bar(x + width/2, best_vals, width,
           color=[RECEPTOR_COLORS[r] for r in receptors], edgecolor="none", label="Best candidate")

    for i, v in enumerate(best_vals):
        improvement = (known_vals[i] - v) / known_vals[i] * 100
        ax.annotate(f"{v:.0f}\n({improvement:+.1f}%)", (i + width/2, v + 15),
                    ha="center", va="bottom", fontsize=5.5)

    ax.set_xticks(x)
    ax.set_xticklabels(receptors)
    ax.set_ylabel("RAPiDock ref2015 score")
    ax.legend(loc="upper left", fontsize=6)
    ax.set_title("a", loc="left", fontweight="bold", fontsize=9)
    ax.set_ylim(0, max(known_vals + best_vals) * 1.18)

    # ---------- Panel b: all candidates distribution ----------
    ax2 = fig.add_subplot(gs[0, 1])
    all_dfs = []
    for rec in receptors:
        df = load_final(rec).dropna(subset=["pos_score"])
        all_dfs.append(df.assign(receptor=rec))
    alld = pd.concat(all_dfs)
    # normalized z-score within receptor (scores not comparable across receptors)
    for rec in receptors:
        mask = alld["receptor"] == rec
        mean = alld.loc[mask, "pos_score"].mean()
        std = alld.loc[mask, "pos_score"].std()
        alld.loc[mask, "pos_norm"] = (alld.loc[mask, "pos_score"] - mean) / std
        alld.loc[mask, "known_norm"] = (KNOWN_POS[rec] - mean) / std

    positions = np.arange(len(receptors))
    for i, rec in enumerate(receptors):
        data = alld[alld["receptor"] == rec]["pos_norm"]
        jitter = np.random.default_rng(i).normal(0, 0.06, len(data))
        ax2.scatter(np.full(len(data), i) + jitter, data, s=14, alpha=0.6,
                    color=RECEPTOR_COLORS[rec], edgecolor="none", zorder=3)
        bp = ax2.boxplot(data, positions=[i], widths=0.45, showfliers=False,
                         patch_artist=True, medianprops=dict(color="black", linewidth=0.8))
        bp["boxes"][0].set_facecolor("white")
        bp["boxes"][0].set_edgecolor(RECEPTOR_COLORS[rec])
        bp["boxes"][0].set_linewidth(0.8)
        # known line
        known_norm = alld[alld["receptor"] == rec]["known_norm"].iloc[0]
        ax2.plot([i - 0.4, i + 0.4], [known_norm, known_norm], color=GRAY,
                 linestyle="--", linewidth=1.0, zorder=4)
    ax2.set_xticks(positions)
    ax2.set_xticklabels(receptors)
    ax2.set_ylabel("ref2015 score (z-scored within receptor)")
    ax2.plot([], [], color=GRAY, linestyle="--", label="Known ligand")
    # 顶部留白, 使图例避开 NK1R 最高离群点 (z≈3.1)
    _v = np.concatenate([alld.loc[alld["receptor"] == r, "pos_norm"].values for r in receptors])
    ax2.set_ylim(_v.min() - 0.25, _v.max() + 0.90)
    ax2.legend(loc="upper right", fontsize=6)
    ax2.set_title("b", loc="left", fontweight="bold", fontsize=9)

    # ---------- Panel c: ΔPos vs ΔNeg quadrant (hero) ----------
    ax3 = fig.add_subplot(gs[1, 0])
    for rec in receptors:
        df = load_final(rec).dropna(subset=["ΔPos", "ΔNeg"])
        ax3.scatter(df["ΔPos"], df["ΔNeg"], s=16, alpha=0.65,
                    color=RECEPTOR_COLORS[rec], label=rec, edgecolor="white", linewidth=0.3)
    ax3.axhline(0, color=GRAY, linewidth=0.6, linestyle="--")
    ax3.axvline(0, color=GRAY, linewidth=0.6, linestyle="--")
    ax3.axvspan(-2500, 0, color=GREEN, alpha=0.08)
    ax3.axhspan(0, 4000, color=GREEN, alpha=0.08)
    # 三个"每受体最优双优"点都挤在原点附近, 逐点标注必然重叠
    # -> 改用金环标记 + 右上空白区彩色列表 (不画引线, 避免穿过数据)
    # 选择标准与全文一致: both = (ΔPos<0 且 frac_safer>0.5), 取 ΔPos 最小者
    best_dual = {}
    for rec in receptors:
        df = load_final(rec).dropna(subset=["ΔPos", "ΔNeg"])
        cand = df[df["both"] == True] if "both" in df else df.iloc[0:0]
        if len(cand) == 0:
            continue
        sel = cand.sort_values("ΔPos").iloc[0]
        best_dual[rec] = sel
        ax3.scatter([sel["ΔPos"]], [sel["ΔNeg"]], s=85, facecolors="none",
                    edgecolors=GOLD, linewidths=1.3, zorder=7)

    ax3.set_xlabel("ΔPos (candidate − known)\nlower = stronger target binding")
    ax3.set_ylabel("ΔNeg (candidate − known)\nhigher = weaker off-target binding")
    ax3.set_title("c", loc="left", fontweight="bold", fontsize=9)
    ax3.set_xlim(-350, 2500)
    ax3.set_ylim(-600, 1800)

    # 右上空白区 (y>850 无任何数据点): 先说明理想象限, 再列出三个最优候选
    ax3.text(0.985, 0.985, "shaded = ideal region:", transform=ax3.transAxes,
             ha="right", va="top", fontsize=5.5, color=DGREEN, fontweight="bold")
    ax3.text(0.985, 0.925, "stronger binding + safer off-target", transform=ax3.transAxes,
             ha="right", va="top", fontsize=5.2, color=DGREEN)
    ax3.text(0.985, 0.845, "best dual-optimal per receptor", transform=ax3.transAxes,
             ha="right", va="top", fontsize=5.5, color=DGRAY, fontweight="bold")
    for k, rec in enumerate(receptors):
        ax3.text(0.985, 0.775 - k * 0.068, f'{rec}   {best_dual[rec]["code"]}',
                 transform=ax3.transAxes, ha="right", va="top", fontsize=5.5,
                 color=RECEPTOR_COLORS[rec])

    # 图例加白底, 防止与右侧散点视觉粘连
    leg3 = ax3.legend(loc="lower right", fontsize=6)
    leg3.get_frame().set_facecolor("white")
    leg3.get_frame().set_alpha(0.85)
    leg3.get_frame().set_linewidth(0)

    # ---------- Panel d: binding modes (NK2R) ----------
    ax4 = fig.add_subplot(gs[1, 1])
    draw_binding_mode_panel(ax4, "NK2R")
    ax4.set_title("d", loc="left", fontweight="bold", fontsize=9)

    fig.savefig(FIGDIR / "Fig2_docking_performance.svg", bbox_inches="tight")
    fig.savefig(FIGDIR / "Fig2_docking_performance.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / "Fig2_docking_performance.tiff", dpi=600, bbox_inches="tight")
    print("[OK] Fig2 saved (4 panels)")

def draw_binding_mode_panel(ax, receptor):
    """Draw known vs best candidate binding pose overlay for a receptor (2D projection)."""
    from Bio.PDB import PDBParser
    pos_dir = RESULTS_DIR / f"{receptor}_peptide_len10"
    known_dir = {"NK1R": "Substance_P", "NK2R": "pep_known", "NK3R": "Neurokinin_B"}[receptor]
    best_dir = {"NK1R": "pep18", "NK2R": "pep27", "NK3R": "pep21"}[receptor]

    prot_pdb = pos_dir / known_dir / f"{known_dir}_protein_raw.pdb"
    if not prot_pdb.exists():
        prot_pdb = pos_dir / known_dir / "pep_known_protein_raw.pdb"
    known_pdb = pos_dir / known_dir / "rank1_ref2015.pdb"
    best_pdb = pos_dir / best_dir / "rank1_ref2015.pdb"

    def ca_trace(pdb_file):
        coords = []
        resnames = []
        for line in open(pdb_file):
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                resnames.append(line[17:20].strip())
        return np.array(coords), resnames

    kc, kr = ca_trace(known_pdb)
    bc, br = ca_trace(best_pdb)

    # project to 2D (use first two PCs for best separation)
    all_c = np.vstack([kc, bc])
    centered = all_c - all_c.mean(axis=0)
    # simple projection on YZ
    def proj(c):
        return c[:, 1], c[:, 2]

    kx, ky = proj(kc)
    bx, by = proj(bc)

    ax.plot(kx, ky, "-o", color=GRAY, linewidth=1.3, markersize=3.5,
            markerfacecolor="white", markeredgewidth=0.7, label=f"Known ({KNOWN_NAMES[receptor]})")
    ax.plot(bx, by, "-o", color=TEAL, linewidth=1.3, markersize=3.5,
            markerfacecolor="white", markeredgewidth=0.7, label=f"{BEST_CODE[receptor]}")
    # annotate termini — 白底避免压在轨迹线上看不清
    _bb = dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.8)
    ax.annotate("N", (kx[0], ky[0]), fontsize=5, color=GRAY, xytext=(0, -9),
                textcoords="offset points", ha="center", zorder=8, bbox=_bb)
    ax.annotate("C", (kx[-1], ky[-1]), fontsize=5, color=GRAY, xytext=(0, -9),
                textcoords="offset points", ha="center", zorder=8, bbox=_bb)
    ax.annotate("N", (bx[0], by[0]), fontsize=5, color=TEAL, xytext=(0, -9),
                textcoords="offset points", ha="center", zorder=8, bbox=_bb)
    ax.annotate("C", (bx[-1], by[-1]), fontsize=5, color=TEAL, xytext=(0, -9),
                textcoords="offset points", ha="center", zorder=8, bbox=_bb)
    ax.set_xticks([])
    ax.set_yticks([])
    # 底部留白, 让 N 端标注不贴到 x 轴说明文字
    _ylo, _yhi = ax.get_ylim()
    ax.set_ylim(_ylo - (_yhi - _ylo) * 0.07, _yhi)
    ax.set_xlabel(f"{receptor} pocket — Cα traces")
    # 图例移到坐标轴外下方, 不再压在轨迹线上
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2,
              fontsize=6, columnspacing=1.2, handletextpad=0.5)

if __name__ == "__main__":
    fig2()
