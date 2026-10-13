#!/usr/bin/env python3
"""
Figure 4 (Main): 三个受体候选肽的序列特征
- 位置特异性氨基酸频率热图 (position × amino acid)
- 显示 PepMLM 学到的受体特异基序:
  NK1R: S/K 开头, 富含 Q/P, C-F 保守
  NK2R: R/K 开头, CC 双半胱氨酸基序
  NK3R: D/E 开头, P-Y 保守, 富含 Q/S/T

Figure contract:
  Core conclusion: 生成器针对不同受体口袋产生不同的氨基酸偏好 (受体特异设计基序)
  Evidence: 32 候选肽/受体的位置-氨基酸频率
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
# ---------------------------------------------------------------------------

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Liberation Sans', 'DejaVu Sans']
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 7
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False


BLUE = "#0F4D92"
TEAL = "#42949E"
VIOLET = "#9A4D8E"

RECEPTORS = {
    "NK1R": {"color": BLUE, "known": "RPKPQQFFGLM"},
    "NK2R": {"color": TEAL, "known": "HKTDSFVGLM"},
    "NK3R": {"color": VIOLET, "known": "DMHDFFVGLM"},
}

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"

def load_peptides(receptor):
    df = pd.read_csv(ANALYSIS / f"{receptor}_summary.csv")
    df = df[~df["peptide_id"].isin(["Substance_P", "pep_known", "Neurokinin_B"])]
    seqs = [s.strip() for s in df["sequence"].tolist() if isinstance(s, str)]
    # keep only 10-mers for position matrix (NK2R has one 11-mer)
    seqs = [s for s in seqs if len(s) == 10]
    return seqs

def compute_freq(receptor):
    seqs = load_peptides(receptor)
    L = len(seqs[0])
    counts = np.zeros((len(AA_ORDER), L))
    for s in seqs:
        for j, aa in enumerate(s):
            if aa in AA_ORDER:
                counts[AA_ORDER.index(aa), j] += 1
    freq = counts / counts.sum(axis=0, keepdims=True)
    return freq, seqs

def fig4():
    fig = plt.figure(figsize=(183/25.4, 110/25.4))
    gs = fig.add_gridspec(2, 3, wspace=0.15, hspace=0.45, height_ratios=[1.4, 1])

    motifs = {}
    for i, (rec, info) in enumerate(RECEPTORS.items()):
        ax = fig.add_subplot(gs[0, i])
        freq, seqs = compute_freq(rec)
        L = freq.shape[1]
        im = ax.imshow(freq, cmap="Blues", aspect="auto", vmin=0, vmax=0.5)
        ax.set_xticks(range(L))
        ax.set_xticklabels([f"{j+1}" for j in range(L)])
        if i == 0:
            ax.set_yticks(range(len(AA_ORDER)))
            ax.set_yticklabels(list(AA_ORDER))
        else:
            ax.set_yticks([])
        ax.set_title(f"{'abc'[i]} — {rec} candidates (n={len(seqs)})", loc="left",
                     fontweight="bold", fontsize=8)
        ax.tick_params(length=2)

        # conserved positions (>40%)
        for j in range(L):
            if freq[:, j].max() > 0.4:
                ax.add_patch(plt.Rectangle((j - 0.5, -0.5), 1, len(AA_ORDER),
                                           fill=False, edgecolor=info["color"], linewidth=1.2))
        # consensus motif
        motif = "".join(AA_ORDER[np.argmax(freq[:, j])] if freq[:, j].max() > 0.3 else "x"
                        for j in range(L))
        motifs[rec] = motif
        ax.set_xlabel(f"consensus: {motif}", fontsize=5.5, color=info["color"])

    cbar_ax = fig.add_axes([0.93, 0.56, 0.012, 0.3])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("Frequency", fontsize=6)
    cbar.ax.tick_params(labelsize=5.5)

    # ---------- Panel d: receptor-specific motif comparison ----------
    axd = fig.add_subplot(gs[1, :])
    # Compare top-2 amino acids per position across receptors
    positions = range(1, 11)
    recs = list(RECEPTORS.keys())
    freqs = {rec: compute_freq(rec)[0] for rec in recs}
    bar_w = 0.25
    x = np.arange(10)
    for i, rec in enumerate(recs):
        top = [AA_ORDER[np.argmax(freqs[rec][:, j])] for j in range(10)]
        topf = [freqs[rec][:, j].max() for j in range(10)]
        axd.bar(x + (i - 1) * bar_w, topf, bar_w, color=RECEPTORS[rec]["color"],
                edgecolor="none", label=rec)
        for j, (aa, f) in enumerate(zip(top, topf)):
            axd.annotate(aa, (x[j] + (i - 1) * bar_w, f + 0.01), ha="center",
                         fontsize=5.5, fontweight="bold", color=RECEPTORS[rec]["color"])
    axd.set_xticks(x)
    axd.set_xticklabels([f"{p}" for p in positions])
    axd.set_xlabel("Position in 10-mer peptide")
    axd.set_ylabel("Frequency of dominant residue")
    axd.set_ylim(0, 0.75)
    axd.legend(loc="upper right", fontsize=6)
    axd.set_title("d — dominant residue per position (receptor-specific motifs)",
                  loc="left", fontweight="bold", fontsize=8)

    fig.savefig(FIGDIR / "Fig4_sequence_features.svg", bbox_inches="tight")
    fig.savefig(FIGDIR / "Fig4_sequence_features.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / "Fig4_sequence_features.tiff", dpi=600, bbox_inches="tight")
    print("[OK] Fig4 saved (4 panels)")

    # Print conserved motifs for manuscript text
    for rec, motif in motifs.items():
        print(f"  {rec} consensus motif: {motif}")

if __name__ == "__main__":
    fig4()
