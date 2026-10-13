#!/usr/bin/env python3
"""
figures4papers house style (Chen Liu, https://github.com/ChenLiu-1996/figures4papers)

从仓库提取的规范:
  PALETTE / DEFAULT_COLORS      — references/api.md
  FigureStyle                   — font_size, axes_linewidth
  标定规则                       — demo 中 font.size 与 linewidth 随画布尺寸缩放:
       大画布 (28 in 宽): font 24, linewidth 3   → 0.857 pt/in, 0.107 lw/in
       紧凑画布 (9–14 in): font 15–16, linewidth 2 → 1.6 pt/in,   0.18 lw/in

应用到本项目 (稿件图幅 183 mm = 7.2 in):
  采用大画布档的密度标定 (font ≈ 6.2 pt, lw ≈ 0.77), 与既有图一致;
  但按 house style 加粗坐标轴与刻度、统一语义配色、并去除残余杂物。

用法:
    from fig_style import PALETTE, apply_style, save_pub, legend_panel
    apply_style()                     # 设定 rcParams
    save_pub(fig, "figures/FigX")     # 导出 svg + pdf + tiff(600 dpi)
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

# ----------------------------------------------------------------------------
# 语义配色 (references/api.md 的 PALETTE)
# ----------------------------------------------------------------------------
PALETTE = {
    "blue_main":      "#0F4D92",
    "blue_secondary": "#3775BA",
    "green_1": "#DDF3DE", "green_2": "#AADCA9", "green_3": "#8BCF8B",
    "red_1":   "#F6CFCB", "red_2":   "#E9A6A1", "red_strong": "#B64342",
    "neutral":   "#CFCECE",
    "highlight": "#FFD700",
    "teal":   "#42949E",
    "violet": "#9A4D8E",
}

# 未显式指定颜色时的默认次序
DEFAULT_COLORS = [PALETTE["blue_main"], PALETTE["green_3"], PALETTE["red_strong"],
                  PALETTE["teal"], PALETTE["violet"], PALETTE["neutral"]]

# 受体配色: 取 house palette 中语义最近的三色 (与稿件 Fig 1–4 的蓝/青/紫一致)
RECEPTOR_COLORS = {
    "NK1R": PALETTE["blue_main"],    # #0F4D92
    "NK2R": PALETTE["teal"],         # #42949E
    "NK3R": PALETTE["violet"],       # #9A4D8E
}

AMBER = "#D9A441"          # 敏感性阈值线 (与稿件 Fig 2b/3 的 dotted 0.75 一致)
GRAY  = PALETTE["neutral"]
DGRAY = "#4D4D4D"

# 稿件图幅下的密度标定 (183 mm = 7.20 in)
_BASE_FONT = 7.0
_BASE_LW   = 1.1           # house style 比默认更粗


def apply_style(base_font=_BASE_FONT, base_lw=_BASE_LW):
    """设定 figures4papers house style 的 rcParams。"""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": base_font,
        "axes.linewidth": base_lw,          # house style: 加粗坐标轴
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.labelsize": base_font,
        "axes.titlesize": base_font + 1,
        "xtick.labelsize": base_font - 1,
        "ytick.labelsize": base_font - 1,
        "xtick.major.width": base_lw * 0.9,
        "ytick.major.width": base_lw * 0.9,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "legend.frameon": False,            # house style: 无边框图例
        "legend.fontsize": base_font - 1,
        "lines.linewidth": base_lw * 1.1,
        "lines.markersize": 4,
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "svg.fonttype": "none",             # SVG 文本可编辑
        "pdf.fonttype": 42,                 # PDF 嵌入 TrueType, 文本可编辑
        "axes.grid": False,                 # house style: 无网格
    })


def save_pub(fig, stem, dpi=600, outdir=None):
    """按 house style 导出: SVG + PDF + 600 dpi TIFF (LZW)。"""
    stem = Path(stem)
    if outdir is not None:
        stem = Path(outdir) / stem.name
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{stem}.svg", bbox_inches="tight")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.tiff", dpi=dpi, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"})
    return stem


def legend_panel(fig, gs_slot, handles, labels, ncol=1, **kw):
    """
    Dedicated legend panel (common-patterns.md):
    在独立子图上 set_axis_off() 后绘制图例, 避免图例遮挡数据。

    fig, gs 由调用方提供; gs_slot 为 (row, col) 或 slice。
    """
    ax = fig.add_subplot(gs_slot)
    ax.set_axis_off()
    ax.legend(handles, labels, ncol=ncol, loc="center", frameon=False, **kw)
    return ax


def sig_bracket(ax, x1, x2, y, p, h=None, fs=5.5, color=None, label=None):
    """
    在两组之间绘制显著性括号并标注 p 值 (Mann-Whitney U)。

    ax   : 目标坐标轴
    x1,x2: 两组在 x 轴上的位置
    y    : 括号横线的数据坐标高度 (通常取两组最大值之上)
    p    : p 值; 也接受预格式化的字符串
    h    : 括号竖直下垂高度 (默认取 y 轴的 3%)
    """
    color = color or DGRAY
    if h is None:
        lo, hi = ax.get_ylim()
        h = (hi - lo) * 0.030
    ax.plot([x1, x1, x2, x2], [y - h, y, y, y - h],
            color=color, linewidth=0.9, solid_capstyle="butt", clip_on=False, zorder=6)
    if isinstance(p, str):
        txt = p
    elif p < 1e-4:
        txt = "P < 0.0001"
    else:
        txt = f"P = {p:.4f}" if p >= 1e-3 else f"P = {p:.2e}"
    if label:
        txt = f"{label}  {txt}"
    ax.annotate(txt, xy=((x1 + x2) / 2, y), xytext=(0, 1.6),
                textcoords="offset points", ha="center", va="bottom",
                fontsize=fs, color=color, zorder=6)
    return txt


def mwu_p(a, b):
    """双侧 Mann-Whitney U 检验 (不假设正态, 适合小样本 docking score 比较)。"""
    from scipy import stats
    a = [v for v in a if v == v]
    b = [v for v in b if v == v]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
