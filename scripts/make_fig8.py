#!/usr/bin/env python3
"""
Figure 8 (Main): ESM-2 潜在空间中的 NK2R 结合肽设计景观
中心: UMAP (拆为 岛屿区 + 主体区 两个并列面板), 四周: 4 个 PDB 结合姿态

布局:
  ┌────────┬───────────────────────┬────────┐
  │ PDB a  │ e1 岛屿 │ e2 主体       │ PDB c  │
  ├────────┤                       ├────────┤
  │ PDB b  │                       │ PDB d  │
  ├────────┴───────────────────────┴────────┤
  │   f 长度    │   g 新颖性   │  h 组成     │
  └─────────────────────────────────────────┘
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from pdb_render import render, parse_atoms
from fig_style import (PALETTE, RECEPTOR_COLORS, AMBER, GRAY, DGRAY as _DG,
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

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Liberation Sans', 'DejaVu Sans']
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 7
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['legend.frameon'] = False

RESULTS = RESULTS_DIR

ORANGE = "#E28E2C"
GOLD = "#C9971C"
BLUE = "#0F4D92"
TEAL = "#3E8E8E"
PURPLE = "#8B5FA8"
GRAY = "#B8B8B8"
DGRAY = "#4D4D4D"
LGRAY = "#E0E0E0"
ISLAND_X = -10.0   # 岛屿/主体分界


# ---------------- PDB 姿态 ----------------
def pose_files(rec_or_dir, pid=None):
    """返回 (protein_raw, peptide_rank1)"""
    if pid is None:                      # 直接给目录
        d = Path(rec_or_dir)
    else:
        d = RESULTS / f"{rec_or_dir}_peptide_len10" / pid
    pr = sorted(d.glob("*_protein_raw.pdb"))
    pe = d / "rank1_ref2015.pdb"
    if pr and pe.exists():
        return pr[0], pe
    return None, None


def pocket_center_of(pep_pdb):
    atoms = parse_atoms(pep_pdb)
    return np.array([a["xyz"] for a in atoms]).mean(axis=0)


def add_pose(ax, pr, pe, title, sub, view_center, color=DGRAY):
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
    try:
        render(ax, pr, pe, title="", pocket_view=True, zoom_radius=19.0,
               view_center=view_center, elev=16, azim=40)
    except Exception as e:
        ax.annotate(f"render failed:\n{e}", xy=(0.5, 0.5), xycoords="axes fraction",
                    ha="center", fontsize=5)
    ax.set_title(title, loc="left", fontweight="bold", fontsize=7, pad=2, color=color)
    if sub:
        ax.annotate(sub, xy=(0.02, 0.01), xycoords="axes fraction",
                    fontsize=5, color=DGRAY, va="bottom", ha="left")


# ---------------- UMAP 绘制 ----------------
def draw_umap(ax, d, pt=4, top_s=34, ours_s=26):
    ax.scatter(d["umap1"], d["umap2"], s=pt, c=LGRAY, alpha=0.6, linewidths=0, zorder=1)
    t10 = d[d["is_top10"] == True]
    ax.scatter(t10["umap1"], t10["umap2"], s=top_s, c=ORANGE,
               edgecolors="white", linewidths=0.6, zorder=4)
    for rec, col, mk in [("ours_NK2R", BLUE, "o"), ("ours_NK1R", TEAL, "^"),
                         ("ours_NK3R", PURPLE, "s")]:
        sub = d[d["group"] == rec]
        if len(sub):
            ax.scatter(sub["umap1"], sub["umap2"], s=ours_s, c=col, marker=mk,
                       edgecolors="white", linewidths=0.4, zorder=5)
    if "is_dual" in d:
        dl = d[d["is_dual"] == True]
        if len(dl):
            ax.scatter(dl["umap1"], dl["umap2"], s=ours_s * 3.2, facecolors="none",
                       edgecolors=GOLD, linewidths=1.2, zorder=6)
    if "in_leaderboard" in d:
        ol = d[d["in_leaderboard"] == True]
        if len(ol):
            ax.scatter(ol["umap1"], ol["umap2"], s=ours_s * 2.4, facecolors="none",
                       edgecolors=ORANGE, linewidths=1.0, zorder=6)


def main():
    df = pd.read_csv(ANALYSIS / "umap_coords.csv")
    main_d = df[df["umap1"] > ISLAND_X].copy()
    isl_d = df[df["umap1"] <= ISLAND_X].copy()

    # ---- 姿态结构文件 (统一相机: 以对照肽口袋质心为参照) ----
    known_pr, known_pe = pose_files("NK2R", "pep_known")
    comp01_pr, comp01_pe = pose_files(RESULTS / "COMP_TOP10_NK2R" / "comp01")
    comp02_pr, comp02_pe = pose_files(RESULTS / "COMP_TOP10_NK2R" / "comp02")
    ours_pr,   ours_pe   = pose_files("NK2R", "pep13")
    vc = pocket_center_of(known_pe)

    # ---- 4 个结合姿态面板对应的 UMAP 坐标 (用于引线定位) ----
    def coord(lbl):
        r = df[df["label"] == lbl].iloc[0]
        return r["umap1"], r["umap2"]

    P = {
        "a": dict(lbl="LB0002",  ax=None, xy=coord("LB0002")),   # 岛屿代表 COMP-02
        "b": dict(lbl="NK2R-13", ax=None, xy=coord("NK2R-13")),  # 本研究双优候选
        "c": dict(lbl="NK2R-K",  ax=None, xy=coord("NK2R-K")),   # 受体匹配对照
        "d": dict(lbl="LB0001",  ax=None, xy=coord("LB0001")),   # 竞赛冠军 COMP-01
    }

    fig = plt.figure(figsize=(183/25.4, 196/25.4))
    # 行: [大 e] / [4 个姿态] / [f g h]
    # 外层 3 段 (上块/图例行/fgh行) 紧凑; 上块内部 (e图与姿态行) 另用大间距,
    # 以避开 e 图的 x 轴刻度标签 —— 避免统一 hspace 把姿态行与图例也拉远
    outer = fig.add_gridspec(3, 1, height_ratios=[3.21, 0.22, 0.62], hspace=0.10)
    top = outer[0].subgridspec(2, 4, height_ratios=[1.75, 0.92],
                               hspace=0.62, wspace=0.30)
    bot = outer[2].subgridspec(1, 4, wspace=0.30)

    # ============ 大 e: 单一 UMAP (断轴) ============
    # 岛屿 (umap1 ~ -21.7..-20.6) 与主体 (2.1..9.4) 在 UMAP-1 上相隔 ~22 单位.
    # 直接同轴绘制会让岛屿缩成不可读的窄条, 故用"断轴"把两段并置为一张图:
    #   共享 y 轴 (可直接比较 UMAP-2), x 轴在空白区间断开并标注断口符号.
    inner = top[0, :].subgridspec(1, 2, width_ratios=[0.80, 3.05], wspace=0.030)
    ax_i = fig.add_subplot(inner[0, 0])   # 岛屿段
    ax_m = fig.add_subplot(inner[0, 1])   # 主体段

    # --- 岛屿段 ---
    draw_umap(ax_i, isl_d, pt=22, top_s=60, ours_s=44)
    ax_i.set_xlim(isl_d["umap1"].min() - 0.10, isl_d["umap1"].max() + 0.28)
    ax_i.set_ylim(main_d["umap2"].min() - 0.35, main_d["umap2"].max() + 0.35)
    ax_i.set_title("e", loc="left", fontweight="bold", fontsize=9, pad=6)
    ax_i.set_xticks([-21.5, -21.0, -20.5])
    ax_i.tick_params(labelsize=5, length=2)
    ax_i.set_ylabel("UMAP-2", fontsize=6)
    ax_i.annotate(f"n = {len(isl_d)}\nlength 7–10 aa", xy=(0.03, 0.97),
                  xycoords="axes fraction", ha="left", va="top", fontsize=5, color=DGRAY)

    # --- 主体段 (共享 y, 隐藏左脊与 y 刻度) ---
    draw_umap(ax_m, main_d, pt=4, top_s=34, ours_s=26)
    ax_m.set_xlim(main_d["umap1"].min() - 0.25, main_d["umap1"].max() + 0.30)
    ax_m.set_ylim(ax_i.get_ylim())
    ax_m.spines["left"].set_visible(False)
    ax_m.set_yticks([])
    ax_m.tick_params(labelsize=5, length=2)
    ax_m.annotate(f"n = {len(main_d)}\nlength 3–40 aa", xy=(0.03, 0.97),
                  xycoords="axes fraction", ha="left", va="top", fontsize=5, color=DGRAY)

    # --- 断口符号 (斜杠) + 共用 x 轴标签 ---
    _d, _kw = 0.012, dict(transform=ax_i.transAxes, color=DGRAY, clip_on=False,
                          linewidth=0.9)
    ax_i.plot((1 - _d, 1 + _d), (-_d * 2.2, +_d * 2.2), **_kw)      # 左段右端断口
    _kw2 = dict(transform=ax_m.transAxes, color=DGRAY, clip_on=False, linewidth=0.9)
    ax_m.plot((-_d, +_d), (-_d * 2.2, +_d * 2.2), **_kw2)           # 右段左端断口
    ax_i.set_xlabel("UMAP-1", fontsize=6, labelpad=1)
    ax_m.set_xlabel("UMAP-1", fontsize=6, labelpad=1)
    # 断口说明贴在断口上方
    ax_m.annotate("axis broken", xy=(0.0, -0.012), xycoords="axes fraction",
                  xytext=(2, -9), textcoords="offset points",
                  fontsize=5, color=DGRAY, ha="left", va="top")
    # 两区小标题贴在各自面板顶部
    ax_i.annotate("Island", xy=(0.5, 1.015), xycoords="axes fraction",
                  ha="center", va="bottom", fontsize=6.5, color=DGRAY, style="italic")
    # Island 简要说明: 置于面板上方空白区 (该区无数据点), 避免压到下方姿态面板
    ax_i.annotate("short-peptide designs (7–10 aa)\nwith non-canonical motifs,\n"
                  "separated from the main space;\n4 of the competition top-10 fall here",
                  xy=(0.0, 0.74), xycoords="axes fraction",
                  ha="left", va="top", fontsize=4.7, color=DGRAY,
                  linespacing=1.45, clip_on=False, zorder=9,
                  bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.85))
    ax_m.annotate("Main sequence space", xy=(0.5, 1.015), xycoords="axes fraction",
                  ha="center", va="bottom", fontsize=6, color=DGRAY, style="italic")

    # 在 e 中给 4 个姿态点加实心端点标记 (识别靠引线, 不靠圈)
    for k, d in P.items():
        xx, yy = d["xy"]
        ax_t = ax_i if xx <= ISLAND_X else ax_m
        ax_t.scatter([xx], [yy], s=16, c="black", zorder=8,
                     edgecolors="white", linewidths=0.6)

    leg = [
        Line2D([], [], marker="o", ls="", ms=3, mfc=LGRAY, mec=LGRAY, label="Leaderboard top-1,000"),
        Line2D([], [], marker="o", ls="", ms=4.5, mfc=ORANGE, mec="white", label="Competition top-10"),
        Line2D([], [], marker="o", ls="", ms=4.5, mfc=BLUE, mec="white", label="This study — NK2R"),
        Line2D([], [], marker="^", ls="", ms=4, mfc=TEAL, mec="white", label="This study — NK1R"),
        Line2D([], [], marker="s", ls="", ms=4, mfc=PURPLE, mec="white", label="This study — NK3R"),
        Line2D([], [], marker="o", ls="", ms=5.5, mfc="none", mec=GOLD, mew=1.2, label="Dual-optimal"),
        Line2D([], [], marker="o", ls="", ms=5, mfc="none", mec=ORANGE, mew=1.0, label="Also in leaderboard"),
    ]
    ax_m.legend(handles=leg, loc="lower right", fontsize=4.8, ncol=2,
                handletextpad=0.3, columnspacing=0.7, borderpad=0.35)

    # ============ 姿态面板: 按目标点在 e 中的 x 位置从左到右排列 ============
    pose_specs = [
        ("a", comp02_pr, comp02_pe, "Island representative\nCOMP-02 / GFPCFYFLM",
         "EC$_{50}$ = 7.93 nM", ORANGE),
        ("b", ours_pr,   ours_pe,   "This study — dual-optimal\nNK2R-13 KSNTQVGFGM",
         "score 660 · frac_safer 0.52", BLUE),
        ("c", known_pr,  known_pe,  "Receptor-matched control\nHKTDSFVGLM",
         "baseline · score 771", DGRAY),
        ("d", comp01_pr, comp01_pe, "Competition champion\nCOMP-01 / Y0002",
         "EC$_{50}$ = 2.07 nM", ORANGE),
    ]
    for col, (key, pr, pe, title, sub, colr) in enumerate(pose_specs):
        axp = fig.add_subplot(top[1, col])
        add_pose(axp, pr, pe, title, sub, vc, colr)
        P[key]["ax"] = axp

    # ---- a–d 姿态面板的残基配色图例 (ClustalX 理化性质分组) ----
    axL = fig.add_subplot(outer[1]); axL.set_axis_off()
    _aa_leg = [
        ("#FF9C00", "G P S T", "small"),
        ("#00BB00", "H K R",   "basic"),
        ("#4D94FF", "F W Y",   "aromatic"),
        ("#9C9C9C", "I L M V", "hydrophobic"),
        ("#E63333", "D E",     "acidic"),
        ("#9933FF", "N Q",     "polar"),
        ("#FF66CC", "C",       "cysteine"),
    ]
    _h = [Line2D([], [], marker="s", ls="", ms=5.5, mfc=c, mec="none",
                 label=f"{aa}  ({d})") for c, aa, d in _aa_leg]
    axL.legend(handles=_h, loc="upper center", bbox_to_anchor=(0.5, 0.98), ncol=7,
               fontsize=5.2, handletextpad=0.35, columnspacing=1.5, borderpad=0.2,
               title="Peptide residue colouring in a–d  (physicochemical class;"
                     " grey ribbon = backbone, blue-grey sticks = pocket residues ≤8 Å)",
               title_fontsize=5.2)
    # (ribbon/sticks/统一相机等细节移入图注, 避免在此与图例条目争位置)

    # ---- 引线: 姿态面板顶边 → e 中的目标点 (代替圈标) ----
    from matplotlib.patches import ConnectionPatch
    for key, d in P.items():
        axp = d["ax"]
        xx, yy = d["xy"]
        ax_t = ax_i if xx <= ISLAND_X else ax_m
        fig.add_artist(ConnectionPatch(
            xyA=(0.985, 1.0), coordsA=axp.transAxes,     # 右上角, 避开左对齐标题
            xyB=(xx, yy), coordsB=ax_t.transData,
            color="#8C8C8C", linewidth=0.55, linestyle="-",
            shrinkA=3, shrinkB=4, zorder=1, alpha=0.85))
        # 面板字母放在标题之上 (左上角), 不与引线交叠
        axp.annotate(key, xy=(0.0, 1.0), xycoords="axes fraction",
                     xytext=(0, 15), textcoords="offset points",
                     ha="left", va="bottom", fontsize=8.5, fontweight="bold",
                     color=colr, zorder=10)

    # ============ 底部分析面板 ============
    # f: 长度 (仅主体区, 避免被岛屿拉伸)
    axf = fig.add_subplot(bot[0, 0:2])
    sc = axf.scatter(main_d["umap1"], main_d["umap2"], c=main_d["length"], s=3.5,
                     cmap="viridis", vmin=5, vmax=25, linewidths=0)
    cb = fig.colorbar(sc, ax=axf, orientation="horizontal", fraction=0.055, pad=0.05,
                      aspect=14)
    cb.set_label("Length (aa)", fontsize=4.8, labelpad=1)
    cb.ax.tick_params(labelsize=4.2, length=1, pad=1)
    axf.set_xticks([]); axf.set_yticks([])
    axf.set_title("f  Length", loc="left", fontweight="bold", fontsize=7)
    for s in axf.spines.values():
        s.set_color(LGRAY)

    # g: 新颖性
    axg = fig.add_subplot(bot[0, 2])
    ours = df[df["group"].str.startswith("ours") & df["nn_dist_lb"].notna()].copy()
    is_dual = ours["is_dual"] == True if "is_dual" in ours else pd.Series(False, index=ours.index)
    d1 = ours[is_dual]["nn_dist_lb"]; d2 = ours[~is_dual]["nn_dist_lb"]
    bins = np.linspace(0, 0.30, 16)
    axg.hist([d2, d1], bins=bins, stacked=True, color=[BLUE, GOLD],
             edgecolor="white", linewidth=0.4,
             label=[f"Other candidates (n={len(d2)})", f"Dual-optimal (n={len(d1)})"])
    axg.axvline(ours["nn_dist_lb"].median(), color=DGRAY, ls="--", lw=0.8)
    axg.annotate(f"median {ours['nn_dist_lb'].median():.2f}", xy=(ours["nn_dist_lb"].median(), 0),
                 xytext=(0.02, 0.94), textcoords="axes fraction", fontsize=5,
                 color=DGRAY, va="top")
    axg.set_xlabel("Distance to nearest leaderboard peptide (UMAP)", fontsize=6)
    axg.set_ylabel("Candidates", fontsize=6)
    axg.tick_params(labelsize=5, length=2)
    axg.legend(fontsize=4.8, loc="upper right")
    axg.set_title("g  Novelty", loc="left", fontweight="bold", fontsize=7)

    # h: 组成
    axh = fig.add_subplot(bot[0, 3])
    counts = {
        "LB\nothers": int((df["group"] == "leaderboard").sum()),
        "Own\nteam": int(df.get("is_own_team", pd.Series(False, index=df.index)).sum()),
        "Ours\nNK2R": int((df["group"] == "ours_NK2R").sum()),
        "Ours\nNK1R/3R": int(df["group"].isin(["ours_NK1R", "ours_NK3R"]).sum()),
    }
    cols = [LGRAY, ORANGE, BLUE, TEAL]
    axh.bar(range(4), list(counts.values()), color=cols, edgecolor="white", linewidth=0.5)
    for i, v in enumerate(counts.values()):
        axh.annotate(str(v), (i, v + 15), ha="center", fontsize=5)
    axh.set_xticks(range(4))
    axh.set_xticklabels(list(counts.keys()), fontsize=5.2)
    axh.set_ylabel("Peptides", fontsize=6)
    axh.tick_params(labelsize=5, length=2)
    axh.set_ylim(0, 1080)
    axh.set_title("h  Composition", loc="left", fontweight="bold", fontsize=7)

    fig.savefig(FIGDIR / "Fig8_latent_space.svg", bbox_inches="tight")
    fig.savefig(FIGDIR / "Fig8_latent_space.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / "Fig8_latent_space.tiff", dpi=600, bbox_inches="tight")
    print("[OK] Fig8 saved")


if __name__ == "__main__":
    main()
