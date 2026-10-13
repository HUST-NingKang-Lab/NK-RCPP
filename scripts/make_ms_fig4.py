#!/usr/bin/env python3
"""
Manuscript Fig. 4 | Position-wise amino-acid frequencies in NK-RCPP-generated peptides.

a–c  Frequencies of the 20 standard amino acids at each position of exact 10-residue
     candidates generated for NK1R, NK2R and NK3R.
     - shared frequency scale 0 → 0.60
     - fixed amino-acid order ACDEFGHIKLMNPQRSTVWY
     - consensus character displayed where the most frequent residue reaches >= 0.30
     - NK2R-02 (11 residues) excluded from the positional analysis → NK2R n=31

论文强调: patterns were not tested for enrichment and do not establish
receptor-binding motifs or mechanisms.
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

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Liberation Sans', 'DejaVu Sans']
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.size'] = 7
plt.rcParams['axes.linewidth'] = 0.8


CONTROL = {"NK1R": "Substance_P", "NK2R": "pep_known", "NK3R": "Neurokinin_B"}
AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
VMAX = 0.60                # 稿件声明的共享色标上限
CONSENSUS_MIN = 0.30       # 稿件声明的一致性阈值
RECS = ["NK1R", "NK2R", "NK3R"]


def freq_matrix(rec):
    """返回 (freq[20 × 10], n_used, n_excluded, consensus string)。仅用恰好 10 残基的候选。"""
    df = pd.read_csv(ANALYSIS / f"{rec}_final_ranking.csv")
    df = df[df["peptide_id"] != CONTROL[rec]].copy()
    seqs = [str(s).strip() for s in df["sequence"]]
    ten = [s for s in seqs if len(s) == 10]
    n_excl = len(seqs) - len(ten)
    counts = np.zeros((len(AA_ORDER), 10))
    for s in ten:
        for j, aa in enumerate(s):
            if aa in AA_ORDER:
                counts[AA_ORDER.index(aa), j] += 1
    freq = counts / max(len(ten), 1)
    cons = "".join(AA_ORDER[np.argmax(freq[:, j])] if freq[:, j].max() >= CONSENSUS_MIN else "x"
                   for j in range(10))
    return freq, len(ten), n_excl, cons


def main():
    fig, axes = plt.subplots(1, 3, figsize=(183/25.4, 62/25.4), sharey=True)
    fig.subplots_adjust(wspace=0.10)

    im = None
    for i, rec in enumerate(RECS):
        ax = axes[i]
        freq, n_use, n_excl, cons = freq_matrix(rec)
        im = ax.imshow(freq, cmap="Blues", aspect="auto", vmin=0, vmax=VMAX)

        ax.set_xticks(range(10))
        ax.set_xticklabels(range(1, 11), fontsize=6)
        if i == 0:
            ax.set_yticks(range(len(AA_ORDER)))
            ax.set_yticklabels(list(AA_ORDER), fontsize=6)
        ax.set_xlabel("Peptide position", fontsize=7)
        ax.tick_params(length=2)

        # 两行标题: 第一行受体与样本量, 第二行 consensus 字符串 (避免与标题叠压)
        line1 = f"{'abc'[i]}  {rec} (n={n_use})"
        if n_excl:
            line1 += f"; {n_excl} excluded"
        ax.set_title(f"{line1}\nconsensus ≥{CONSENSUS_MIN:.2f}: {cons}",
                     loc="left", fontweight="bold", fontsize=8, color=C[rec])

    cax = fig.add_axes([0.92, 0.16, 0.014, 0.68])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Position frequency", fontsize=7)
    cb.set_ticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    cb.ax.tick_params(labelsize=6, length=2)
    cb.outline.set_linewidth(0.6)

    fig.savefig(FIGDIR / "Fig4_audit.svg", bbox_inches="tight")
    fig.savefig(FIGDIR / "Fig4_audit.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / "Fig4_audit.tiff", dpi=600, bbox_inches="tight")
    print("[OK] Fig4_audit saved")
    for rec in RECS:
        _, n, ne, cons = freq_matrix(rec)
        print(f"   {rec}: n={n} excluded={ne} consensus={cons}")


if __name__ == "__main__":
    main()
