#!/usr/bin/env python3
"""
Supplementary Figures:
- Fig S1: 每个受体全部候选肽的正向对接分数排名 (vs 已知肽)
- Fig S2: 每个受体全部候选肽的负向配对 ΔNeg 分布 (boxplot)
- Fig S3: 脱靶蛋白的蛋白类型分类 (bar, 3 receptors)
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
# ---------------------------------------------------------------------------
apply_style()

BLUE = PALETTE["blue_main"]
TEAL = PALETTE["teal"]
VIOLET = PALETTE["violet"]
GREEN = "#8BCF8B"
RED = "#B64342"
GOLD = "#E2A93B"


RECEPTORS = {"NK1R": BLUE, "NK2R": TEAL, "NK3R": VIOLET}
KNOWN_POS = {"NK1R": 574.39, "NK2R": 771.13, "NK3R": 2240.54}
KNOWN_NAMES = {"NK1R": "Substance P", "NK2R": "NKA-derived", "NK3R": "Neurokinin B"}

def code_for(rec, pid):
    if pid in ["Substance_P", "pep_known", "Neurokinin_B"]:
        return f"{rec}-K"
    return f"{rec}-{int(pid.replace('pep','')):02d}"

def figS1():
    """Full ranking of all candidates per receptor (by code).

    每受体 33 行, 逐行标注需要 ~2.3mm 行距(5.5pt 字), 故:
      - 图放宽到双栏 (183mm) 并加高, 使行距达到可读水平
      - y 轴隔行标注 (每 2 行一个代号), 避免文字互相叠压
      - 已知肽说明文字放在顶部预留的空白行, 不再压在条形上
    """
    fig, axes = plt.subplots(3, 1, figsize=(183/25.4, 210/25.4), sharex=False)
    fig.subplots_adjust(hspace=0.32)
    for ax, (rec, color) in zip(axes, RECEPTORS.items()):
        df = pd.read_csv(ANALYSIS / f"{rec}_summary.csv")
        df = df.sort_values("pos_score")
        codes = [code_for(rec, pid) for pid in df["peptide_id"]]
        scores = df["pos_score"].tolist()
        ids = df["peptide_id"].tolist()

        colors = []
        for pid, s in zip(ids, scores):
            if pid in ["Substance_P", "pep_known", "Neurokinin_B"]:
                colors.append(GRAY)
            else:
                colors.append(color)
        y = np.arange(len(codes))
        ax.barh(y, scores, color=colors, height=0.78, edgecolor="none")
        ax.axvline(KNOWN_POS[rec], color=GRAY, linestyle="--", linewidth=0.8)

        # 隔行标注: 奇数行留空, 保证相邻标签间距 >= 2 倍行距
        ax.set_yticks(y[::2])
        ax.set_yticklabels([codes[i] for i in range(0, len(codes), 2)], fontsize=5.5)
        for i in range(0, len(codes), 2):          # 留空行画短刻度, 便于对位
            ax.axhline(i + 0.5, color="#EEEEEE", linewidth=0.5, zorder=0)
        ax.set_ylim(len(codes) - 0.5, -2.2)        # 顶部留 2 行放说明文字
        ax.text(KNOWN_POS[rec], -1.6, f"{KNOWN_NAMES[rec]} ({rec}-K)", fontsize=6,
                color=DGRAY, ha="center", va="center", zorder=5,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.9))
        ax.set_xlabel("RAPiDock ref2015 score")
        ax.set_title(f"{rec} — positive docking (lower = better)", loc="left",
                     fontweight="bold", fontsize=8)
        ax.tick_params(length=2, labelsize=6)
    fig.savefig(FIGDIR / "FigS1_positive_ranking.svg", bbox_inches="tight")
    fig.savefig(FIGDIR / "FigS1_positive_ranking.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / "FigS1_positive_ranking.tiff", dpi=600, bbox_inches="tight")
    print("[OK] FigS1 saved")

def figS2():
    """ΔNeg distribution per candidate (per-protein paired deltas vs known).

    同 FigS1: 33 行逐行标注过密 -> 加宽加高 + 隔行标注。
    """
    fig, axes = plt.subplots(3, 1, figsize=(183/25.4, 210/25.4))
    fig.subplots_adjust(hspace=0.32)
    for ax, (rec, color) in zip(axes, RECEPTORS.items()):
        df = pd.read_csv(ANALYSIS / f"{rec}_pairwise_selectivity.csv")
        df = df.sort_values("delta_median", ascending=False)
        df = df.dropna(subset=["delta_mean"])
        # resolve code from sequence via summary table
        summ = pd.read_csv(ANALYSIS / f"{rec}_summary.csv")
        seq2pid = dict(zip(summ["sequence"].str.strip(), summ["peptide_id"]))
        seqs = df["sequence"].tolist()
        codes = [code_for(rec, seq2pid.get(s, "pep0")) for s in seqs]
        medians = df["delta_median"].tolist()
        fracs = df["frac_safer"].tolist()
        y = np.arange(len(codes))
        ax.barh(y, medians, color=[color if f > 0.5 else GRAY for f in fracs],
                height=0.78, edgecolor="none")
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_yticks(y[::2])
        ax.set_yticklabels([codes[i] for i in range(0, len(codes), 2)], fontsize=5.5)
        ax.set_ylim(len(codes) - 0.5, -0.9)
        ax.set_xlabel("ΔNeg median (candidate − known), per off-target protein")
        ax.set_title(f"{rec} — off-target paired deltas (right of 0 = safer)", loc="left",
                     fontweight="bold", fontsize=8)
        ax.tick_params(length=2, labelsize=6)
    fig.savefig(FIGDIR / "FigS2_delta_neg_distribution.svg", bbox_inches="tight")
    fig.savefig(FIGDIR / "FigS2_delta_neg_distribution.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / "FigS2_delta_neg_distribution.tiff", dpi=600, bbox_inches="tight")
    print("[OK] FigS2 saved")

def figS3():
    """Off-target protein class composition."""
    # Classify by known GPCR families from PDB titles
    classification = {
        "NK1R": {"Opioid": 3, "Chemokine": 4, "Biogenic amine": 8, "Peptide GPCR": 9,
                 "Other GPCR": 11, "Glycosyltransferase/synthase": 31},
        "NK2R": {"Chemokine": 2, "Biogenic amine": 4, "Peptide GPCR": 5,
                 "Other GPCR": 7, "Glycosyltransferase/synthase": 9},
        "NK3R": {"Opioid": 3, "Chemokine": 3, "Biogenic amine": 7, "Peptide GPCR": 10,
                 "Other GPCR": 9, "Glycosyltransferase/synthase": 28},
    }
    # Fallback: auto-count from PDB files if available
    try:
        classification = classify_offtargets()
    except Exception as e:
        print(f"[WARN] auto-classification failed ({e}), using manual table")

    fig, axes = plt.subplots(1, 3, figsize=(183/25.4, 60/25.4))
    fig.subplots_adjust(wspace=0.4)
    order = ["Biogenic amine", "Peptide GPCR", "Chemokine", "Opioid", "Other GPCR",
             "Glycosyltransferase/synthase"]
    for ax, (rec, counts) in zip(axes, classification.items()):
        cats = [c for c in order if c in counts]
        vals = [counts[c] for c in cats]
        ax.bar(range(len(cats)), vals, color=RECEPTORS[rec], width=0.6, edgecolor="none")
        ax.set_xticks(range(len(cats)))
        ax.set_xticklabels(cats, rotation=30, ha="right", fontsize=5)
        ax.set_title(f"{rec} off-targets (n={sum(vals)})", loc="left", fontweight="bold", fontsize=7.5)
        ax.set_ylabel("Proteins")
        ax.tick_params(length=2)
    fig.savefig(FIGDIR / "FigS3_offtarget_composition.svg", bbox_inches="tight")
    fig.savefig(FIGDIR / "FigS3_offtarget_composition.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / "FigS3_offtarget_composition.tiff", dpi=600, bbox_inches="tight")
    print("[OK] FigS3 saved")

def classify_offtargets():
    """Auto-classify off-target proteins (only those used in negative docking)
    by keywords from original PDB TITLE lines."""
    import glob, os, pandas as pd
    keyword_map = {
        "Biogenic amine": ["DOPAMINE", "SEROTONIN", "5-HT", "HISTAMINE", "MUSCARINIC",
                           "ACETYLCHOLINE", "ADRENERGIC", "OPSIN", "RHODOPSIN"],
        "Chemokine": ["CHEMOKINE", "CCR", "CXCR"],
        "Opioid": ["OPIOID", "NOCICEPTIN"],
        "Peptide GPCR": ["NEUROPEPTIDE", "OREXIN", "NEUROTENSIN", "GLP-1", "GHRH",
                         "GASTRIN", "MELANOCORTIN", "PROTEASE-ACTIVATED", "NEUROKININ",
                         "BOMBESIN", "NPY", "PAR1", "PAR2", "GLUCAGON", "SOMATOSTATIN",
                         "CORTICOTROPIN", "VASOPRESSIN", "OXYTOCIN", "ANGLOTENSIN",
                         "ENDOTHELIN", "GHRELIN", "KISSPEPTIN", "BRADYKININ", "TACHYKININ"],
        "Other GPCR": ["P2Y", "CB1", "CANNABINOID", "CRTH2", "LEUKOTRIENE", "GPR",
                       "SMOOTHENED", "ADENOSINE", "NIACIN", "OLFACTORY", "LYSOPHOSPHATIDIC",
                       "CALCIUM-SENSING", "EBI2", "GPR183", "HCA2", "MC4R", "MGLU",
                       "SIGNALING PROTEIN", "7TM", "RECEPTOR", "GPCR"],
        "Glycosyltransferase/synthase": ["GLYCOSYLTRANSFERASE", "GLYCOGEN SYNTHASE",
                                         "SUCROSE", "TREHALOSE", "GLUCOSAMINE", "PIMB",
                                         "MNAA", "WECB", "GTFA", "MANNOSYL", "ALDOLASE",
                                         "EPIMERASE", "WBAZ", "WAAG", "WSAF", "TRET",
                                         "MSHA", "STARCH", "BABSHA", "PIMA", "GTFB",
                                         "GLYE", "POLYSIALYLTRANSFERASE", "VLDE", "SPS",
                                         "UDP", "SYNTHASE", "TRANSFERASE"],
    }
    result = {}
    for rec in RECEPTORS:
        # only proteins used in negative docking
        neg_inputs = RAPIDOCK / f"{rec.lower()}_neg_inputs"
        proteins = set()
        f1 = neg_inputs / "peptide01_docking_input.csv"
        if f1.exists():
            df1 = pd.read_csv(f1)
            proteins.update(c.split("_peptide")[0] for c in df1["complex_name"])
        fk = neg_inputs / "known_peptide_docking_input.csv"
        if fk.exists():
            dfk = pd.read_csv(fk)
            proteins.update(c.split("_known")[0] for c in dfk["complex_name"])

        orig_dir = DOCKING / f"{rec.lower()}_foldseek" / "pdb_files"
        counts = {k: 0 for k in keyword_map}
        for pdb_chain in proteins:
            pdb_id = pdb_chain.split("_")[0]
            pdb = orig_dir / f"{pdb_id}.pdb"
            try:
                title = " ".join(l[10:].strip() for l in open(pdb) if l.startswith("TITLE")).upper()
            except Exception:
                title = ""
            assigned = False
            for cat, kws in keyword_map.items():
                if any(kw in title for kw in kws):
                    counts[cat] += 1
                    assigned = True
                    break
            if not assigned:
                counts["Other GPCR"] += 1
        # remove empty categories
        result[rec] = {k: v for k, v in counts.items() if v > 0}
        print(f"  {rec}: {result[rec]}")
    return result

if __name__ == "__main__":
    figS1()
    figS2()
    figS3()
