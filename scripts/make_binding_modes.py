#!/usr/bin/env python3
"""
结合模式结构图生成器
对每个受体: 已知肽 vs 最优候选肽在受体口袋中的对接姿态
- 蛋白: 半透明表面 (Cα 散点), 口袋残基 (距离肽<8Å) 高亮
- 肽: Cα 骨架 trace + 侧链接触
输出: FigS4 (三受体 × 已知/候选 6 panels)
"""
import sys
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).parent))
from fig_style import (PALETTE, RECEPTOR_COLORS, AMBER, GRAY, DGRAY,
                       apply_style, save_pub)
apply_style()
from Bio.PDB import PDBParser, NeighborSearch

sys.path.insert(0, str(Path(__file__).parent))
from pdb_render import AA3 as AA3TO1     # 'ALA' -> 'A'

# --- paths ---------------------------------------------------------------
# Repository-relative: the scripts run from anywhere inside a checkout.
# Raw RAPiDock run directories and the competition workbook are not shipped;
# set NK_RCPP_RAW / NK_RCPP_COMP_XLSX (or NK_RCPP_BASE) to point at them.
# See scripts/README.md.
BASE = Path(os.environ.get("NK_RCPP_BASE", Path(__file__).resolve().parent.parent))
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

BLUE = PALETTE["blue_main"]
TEAL = PALETTE["teal"]
VIOLET = PALETTE["violet"]
DARKGRAY = "#606060"
RED = "#B64342"


# For each receptor: (known_pep_dir, best_pep_dir, color)
COMPARISONS = {
    "NK1R": {"known": "Substance_P", "best": "pep18", "color": BLUE},
    "NK2R": {"known": "pep_known", "best": "pep27", "color": TEAL},
    "NK3R": {"known": "Neurokinin_B", "best": "pep21", "color": VIOLET},
}

def parse_pdb_coords(pdb_file):
    """Parse Cα coordinates: returns list of (resname, chain, resnum, ca_coord)."""
    residues = []
    for line in open(pdb_file):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            resname = line[17:20].strip()
            chain = line[21]
            resnum = int(line[22:26])
            coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            residues.append((resname, chain, resnum, coord))
    return residues

def parse_all_atoms(pdb_file):
    """Parse all atoms for neighbor search."""
    atoms = []
    for line in open(pdb_file):
        if line.startswith("ATOM"):
            name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21]
            resnum = int(line[22:26])
            coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            atoms.append((name, resname, chain, resnum, coord))
    return atoms

def get_peptide_ca(pdb_file):
    """Get peptide Cα trace."""
    return parse_pdb_coords(pdb_file)

def find_pocket_residues(protein_atoms, peptide_atoms, cutoff=8.0):
    """Find protein residues within cutoff of peptide atoms."""
    pep_coords = np.array([a[4] for a in peptide_atoms])
    prot_residues = set()
    for a in protein_atoms:
        d = np.min(np.linalg.norm(pep_coords - a[4], axis=1))
        if d < cutoff:
            prot_residues.add((a[1], a[2], a[3]))  # (resname, chain, resnum)
    return prot_residues

def project_2d(coords, axis="yz"):
    """Project 3D coords to 2D for plotting."""
    if axis == "yz":
        return coords[:, 1], coords[:, 2]
    elif axis == "xy":
        return coords[:, 0], coords[:, 1]
    else:
        return coords[:, 0], coords[:, 2]

def draw_binding_mode(ax, protein_pdb, peptide_pdb, color, title, axis="yz"):
    """Draw a 2D projection binding mode."""
    prot_ca = parse_pdb_coords(protein_pdb)
    prot_atoms = parse_all_atoms(protein_pdb)
    pep_atoms = parse_all_atoms(peptide_pdb)
    pep_ca = parse_pdb_coords(peptide_pdb)

    # pocket residues
    pocket = find_pocket_residues(prot_atoms, pep_atoms, cutoff=8.0)
    pocket_ca = [r for r in prot_ca if (r[0], r[1], r[2]) in pocket]

    prot_coords = np.array([r[3] for r in prot_ca])
    pocket_coords = np.array([r[3] for r in pocket_ca])
    pep_coords = np.array([r[3] for r in pep_ca])

    # protein background (faint)
    px, py = project_2d(prot_coords, axis)
    ax.scatter(px, py, s=3, color=GRAY, alpha=0.25, linewidths=0)

    # pocket residues (medium)
    if len(pocket_coords):
        kx, ky = project_2d(pocket_coords, axis)
        ax.scatter(kx, ky, s=10, color=color, alpha=0.5, linewidths=0)

    # peptide Cα trace
    sx, sy = project_2d(pep_coords, axis)
    ax.plot(sx, sy, "-o", color="black", linewidth=1.2, markersize=4,
            markerfacecolor="white", markeredgewidth=0.8, zorder=5)
    # 残基标签: 三字母太宽会互相叠压 -> 单字母, 并沿局部骨架法线向"背离骨架中心"一侧偏移,
    # 使残基密集处(如近乎垂直的骨架)的标签也互不相撞
    _P = np.column_stack([sx, sy])
    _ctr = _P.mean(axis=0)
    _bb = dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.75)
    for i, r in enumerate(pep_ca):
        j0, j1 = max(0, i - 1), min(len(_P) - 1, i + 1)
        d = _P[j1] - _P[j0]
        n = np.array([-d[1], d[0]])
        nrm = np.linalg.norm(n)
        n = n / nrm if nrm > 1e-6 else np.array([0.0, 1.0])
        if np.dot(n, _P[i] - _ctr) < 0:
            n = -n
        aa1 = AA3TO1.get(r[0], "X")
        ax.annotate(aa1, (sx[i], sy[i]), xytext=(n[0] * 7, n[1] * 7),
                    textcoords="offset points", fontsize=4.5, ha="center", va="center",
                    zorder=7, color=DARKGRAY, bbox=_bb)

    ax.set_title(title, loc="left", fontweight="bold", fontsize=7.5)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

def fig_binding_modes():
    fig, axes = plt.subplots(3, 2, figsize=(183/25.4, 140/25.4))
    fig.subplots_adjust(wspace=0.05, hspace=0.25)

    for i, (rec, info) in enumerate(COMPARISONS.items()):
        pos_dir = RESULTS_DIR / f"{rec}_peptide_len10"
        protein_pdb = pos_dir / info["known"] / f"{info['known']}_protein_raw.pdb"
        if not protein_pdb.exists():
            # NK2R known dir is pep_known
            protein_pdb = pos_dir / info["known"] / "pep_known_protein_raw.pdb"
        known_pdb = pos_dir / info["known"] / "rank1_ref2015.pdb"
        best_pdb = pos_dir / info["best"] / "rank1_ref2015.pdb"

        # Known peptide
        ax1 = axes[i, 0]
        title_known = {"NK1R": "Substance P", "NK2R": "NKA-derived", "NK3R": "Neurokinin B"}[rec]
        draw_binding_mode(ax1, protein_pdb, known_pdb, GRAY, f"{rec} — {title_known}", axis="xy")
        ax1.set_ylabel(rec, fontweight="bold", fontsize=8)

        # Best candidate
        ax2 = axes[i, 1]
        best_seq = {"NK1R": "SPRRQLGFGL", "NK2R": "PSNTQCCFGL", "NK3R": "DEPTQQCIGF"}[rec]
        draw_binding_mode(ax2, protein_pdb, best_pdb, info["color"], f"{rec} — {best_seq}", axis="xy")

    fig.savefig(FIGDIR / "FigS4_binding_modes.svg", bbox_inches="tight")
    fig.savefig(FIGDIR / "FigS4_binding_modes.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / "FigS4_binding_modes.tiff", dpi=600, bbox_inches="tight")
    print("[OK] FigS4 binding modes saved")

if __name__ == "__main__":
    fig_binding_modes()
