#!/usr/bin/env python3
"""
PDB 结合姿态 3D 渲染器 (纯 Bio.PDB + matplotlib, 无需 PyMOL)

风格: 蛋白骨架 ribbon (按深度着色/调透明度) + 口袋残基 sticks + 肽段高亮 tube/sticks
投影: 正交投影 + 画家算法(按深度排序绘制) 实现遮挡
"""
import numpy as np
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle
from pathlib import Path

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

# ---------- 解析 ----------
AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
       'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
       'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
AA1_COLOR = {  # 侧链按理化性质着色 (与 Fig7 一致的 ClustalX 风格)
    'G':'#FF9C00','P':'#FF9C00','S':'#FF9C00','T':'#FF9C00',
    'H':'#00BB00','K':'#00BB00','R':'#00BB00',
    'F':'#4D94FF','W':'#4D94FF','Y':'#4D94FF',
    'I':'#9C9C9C','L':'#9C9C9C','M':'#9C9C9C','V':'#9C9C9C',
    'A':'#C8C8C8','D':'#E63333','E':'#E63333',
    'N':'#9933FF','Q':'#9933FF','C':'#FF66CC',
}


def parse_atoms(pdb):
    """-> list of dict(atoms)"""
    atoms = []
    for line in open(pdb):
        if line.startswith(("ATOM", "HETATM")):
            try:
                atoms.append({
                    "name": line[12:16].strip(),
                    "res": line[17:20].strip(),
                    "chain": line[21],
                    "resid": int(line[22:26]),
                    "xyz": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                    "elem": line[76:78].strip() or line[12:16].strip()[0],
                })
            except ValueError:
                continue
    return atoms


def ca_trace(atoms, chain):
    """Cα 坐标序列 (按残基顺序)"""
    out = []
    for a in atoms:
        if a["chain"] == chain and a["name"] == "CA":
            out.append((a["resid"], a["res"], a["xyz"]))
    return out


def residues_of(atoms, chain):
    d = {}
    for a in atoms:
        if a["chain"] == chain:
            d.setdefault(a["resid"], {"res": a["res"], "atoms": []})["atoms"].append(a)
    return d


def smooth_polyline(pts, n=8):
    """Catmull-Rom 样条平滑, 让骨架更接近 ribbon"""
    if len(pts) < 3:
        return np.array(pts)
    P = np.array(pts, dtype=float)
    ext = np.vstack([P[0], P, P[-1]])
    out = []
    for i in range(1, len(ext) - 2):
        p0, p1, p2, p3 = ext[i-1], ext[i], ext[i+1], ext[i+2]
        for t in np.linspace(0, 1, n, endpoint=False):
            t2, t3 = t*t, t*t*t
            out.append(0.5*((2*p1) + (-p0+p2)*t + (2*p0-5*p1+4*p2-p3)*t2 + (-p0+3*p1-3*p2+p3)*t3))
    out.append(P[-1])
    return np.array(out)


def setup_camera(pts, elev=18, azim=35, center=None):
    """返回投影函数 (正交投影到屏幕坐标)"""
    e, a = np.radians(elev), np.radians(azim)
    fwd = np.array([np.cos(e)*np.cos(a), np.cos(e)*np.sin(a), np.sin(e)])
    right = np.array([-np.sin(a), np.cos(a), 0.0])
    up = np.cross(fwd, right)
    if center is None:
        center = pts.mean(axis=0)
    center = np.asarray(center, dtype=float)
    def project(P):
        Q = np.atleast_2d(P) - center
        return np.column_stack([Q @ right, Q @ up]), (Q @ fwd)
    return project


def render(ax, protein_pdb, peptide_pdb, title="", elev=18, azim=35,
           pocket_cut=8.0, view_center=None, zoom=1.0, pocket_view=True,
           zoom_radius=18.0):
    """渲染一个结合姿态

    pocket_view: True 时裁剪到肽段周围 zoom_radius(Å) 的窗口 (突出结合口袋)
    """
    prot = parse_atoms(protein_pdb)
    pep = parse_atoms(peptide_pdb)

    if pocket_view and pep:
        pep_c = np.array([a["xyz"] for a in pep]).mean(axis=0)
        keep = [a for a in prot if np.linalg.norm(a["xyz"] - pep_c) < zoom_radius]
        prot = keep if len(keep) > 20 else prot

    # ---- 肽段 ----
    pep_chain = pep[0]["chain"] if pep else "A"
    pep_ca = ca_trace(pep, pep_chain)
    pep_res = residues_of(pep, pep_chain)
    pep_coords = np.array([a["xyz"] for a in pep]) if pep else np.zeros((1, 3))

    # ---- 口袋残基 (距肽 < pocket_cut) ----
    pocket_ids = set()
    for a in prot:
        if a["name"] in ("CA", "CB", "CG", "CD", "CZ", "NZ", "OH", "NE2"):
            d = np.linalg.norm(pep_coords - a["xyz"], axis=1).min()
            if d < pocket_cut:
                pocket_ids.add((a["chain"], a["resid"]))
    prot_res = residues_of(prot, prot_chain_of(prot))

    # ---- 相机: 以肽段质心为中心 ----
    if pocket_view and len(pep_coords) > 1:
        center = pep_coords.mean(axis=0)
    else:
        center = np.vstack([np.array([a["xyz"] for a in prot]), pep_coords]).mean(axis=0)
    if view_center is not None:
        center = np.asarray(view_center, dtype=float)
    allpts = np.vstack([np.array([a["xyz"] for a in prot]), pep_coords])
    project = setup_camera(allpts, elev, azim, center=center)

    # ================= 绘制 =================
    # 1) 蛋白骨架: 按链分段平滑, 深度决定透明度/线宽
    chains = sorted(set(a["chain"] for a in prot))
    for ch in chains:
        ca = ca_trace(prot, ch)
        if len(ca) < 4:
            continue
        pts = smooth_polyline([c for _, _, c in ca], n=6)
        xy, depth = project(pts)
        dmin, dmax = depth.min(), depth.max()
        dnorm = (depth - dmin) / max(dmax - dmin, 1e-6)
        # 分段绘制 (每段一个颜色/透明度)
        segs, cols, lws = [], [], []
        for i in range(len(pts) - 1):
            segs.append([xy[i], xy[i+1]])
            # 远的更淡、更细
            alpha = 0.25 + 0.55 * dnorm[i]
            cols.append((*matplotlib.colors.to_rgb("#B9C6D6"), alpha))
            lws.append(1.0 + 2.2 * dnorm[i])
        ax.add_collection(LineCollection(segs, colors=cols, linewidths=lws,
                                         capstyle="round", zorder=1))

    # 2) 口袋残基 sticks (细线)
    segs, cols, lws, deps = [], [], [], []
    for ch, rid in sorted(pocket_ids):
        res = prot_res.get(rid) if ch == prot_chain_of(prot) else None
        if not res:
            continue
        rid_d = {a["resid"]: a for a in res["atoms"]}
        # 找侧链原子 (非 N/CA/C/O)
        side = [a for a in res["atoms"] if a["name"] not in ("N", "CA", "C", "O", "OXT")]
        if not side:
            continue
        base = next((a for a in res["atoms"] if a["name"] == "CA"), None)
        if base is None:
            continue
        chain_pts = [base["xyz"]] + [a["xyz"] for a in sorted(side, key=lambda z: z["name"])]
        if len(chain_pts) < 2:
            continue
        P = np.array(chain_pts)
        xy, depth = project(P)
        dm = depth.mean()
        for i in range(len(P) - 1):
            segs.append([xy[i], xy[i+1]])
            cols.append((*matplotlib.colors.to_rgb("#7A8B99"), 0.60))
            lws.append(1.6)
            deps.append(dm)

    # 3) 肽段: 主链 tube (按残基着色) + 侧链
    pep_pts = [c for _, _, c in pep_ca]
    if len(pep_pts) >= 2:
        sp = smooth_polyline(pep_pts, n=10)
        xy, depth = project(sp)
        # 按最近的 Cα 决定颜色
        ca_xy, _ = project(np.array(pep_pts))
        colors = [AA1_COLOR.get(AA3.get(r, "A"), "#888888") for _, r, _ in pep_ca]
        for i in range(len(sp) - 1):
            # 找最近 Cα 索引
            j = int(np.argmin(np.linalg.norm(ca_xy - xy[i], axis=1)))
            c = colors[min(j, len(colors)-1)]
            ax.plot([xy[i, 0], xy[i+1, 0]], [xy[i, 1], xy[i+1, 1]],
                    color=c, linewidth=5.6, solid_capstyle="round",
                    alpha=0.95, zorder=10 + i * 1e-3)
        # 侧链 sticks
        for rid, resname, ca_xyz in pep_ca:
            res = pep_res.get(rid)
            if not res: continue
            side = [a for a in res["atoms"] if a["name"] not in ("N", "CA", "C", "O", "OXT")]
            if not side: continue
            P = np.array([ca_xyz] + [a["xyz"] for a in sorted(side, key=lambda z: z["name"])])
            xy2, _ = project(P)
            col = AA1_COLOR.get(AA3.get(resname, "A"), "#888888")
            for i in range(len(P) - 1):
                ax.plot([xy2[i,0], xy2[i+1,0]], [xy2[i,1], xy2[i+1,1]],
                        color=col, linewidth=2.0, solid_capstyle="round",
                        alpha=0.9, zorder=20)

    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, loc="left", fontweight="bold", fontsize=8)


def prot_chain_of(atoms):
    """蛋白主链 (残基数最多的链)"""
    from collections import Counter
    c = Counter(a["chain"] for a in atoms if a["name"] == "CA")
    return c.most_common(1)[0][0] if c else "A"


if __name__ == "__main__":
    RESULTS = RESULTS_DIR
    fig, axes = plt.subplots(1, 3, figsize=(180/25.4, 65/25.4))
    cases = [
        (RESULTS/"NK2R_peptide_len10/pep_known/pep_known_protein_raw.pdb",
         RESULTS/"NK2R_peptide_len10/pep_known/rank1_ref2015.pdb", "Known (NKA-derived)"),
        (RESULTS/"NK2R_peptide_len10/pep27/pep27_protein_raw.pdb",
         RESULTS/"NK2R_peptide_len10/pep27/rank1_ref2015.pdb", "This study: NK2R-27"),
        (RESULTS/"COMP_TOP10_NK2R/comp01/comp01_protein_raw.pdb",
         RESULTS/"COMP_TOP10_NK2R/comp01/rank1_ref2015.pdb", "Competition: COMP-01"),
    ]
    for ax, (pp, qp, t) in zip(axes, cases):
        if pp.exists() and qp.exists():
            render(ax, pp, qp, title=t)
        else:
            ax.text(0.5, 0.5, f"missing:\n{pp.name}", ha="center", fontsize=6)
    fig.savefig(str(FIGDIR / "_pdb_render_test.png"),
                dpi=300, bbox_inches="tight", facecolor="white")
    print("[OK] test render saved")
