#!/usr/bin/env python3
"""
Fig 8 步骤1: UMAP 降维
输入: analysis/esm_emb_l2.npy (L2归一化的 ESM2-650M 嵌入), analysis/umap_input.csv
输出: analysis/umap_coords.csv (含 UMAP1/2 坐标)
"""
import numpy as np
import os
import pandas as pd
from pathlib import Path
import umap

# --- paths ---------------------------------------------------------------
# Repository-relative: the scripts run from anywhere inside a checkout.
# Raw RAPiDock run directories and the competition workbook are not shipped;
# set NK_RCPP_RAW / NK_RCPP_COMP_XLSX (or NK_RCPP_BASE) to point at them.
# See scripts/README.md.
BASE = Path(os.environ.get("NK_RCPP_BASE", Path(__file__).resolve().parent.parent))
ANALYSIS = BASE / "step5_analysis" / "output"
# ---------------------------------------------------------------------------


def load_l2_embeddings():
    """
    载入 L2 归一化嵌入 (没有则从原始嵌入生成)。

    必须归一化的原因: 原始 ESM 末层表示的范数跨度达 600 倍 (0.003-1.762),
    低复杂度序列 (如 APAQMLF) 的范数趋近 0。直接用原始向量做 UMAP,
    欧氏距离会被模长主导而非序列语义。详见 UMAP_FIGURE_PLAN.md §3。
    """
    l2_path = ANALYSIS / "esm_emb_l2.npy"
    if l2_path.exists():
        print(f"[INFO] 载入已缓存的 {l2_path.name}")
        return np.load(l2_path)
    raw = np.load(ANALYSIS / "esm_emb.npy")
    norm = np.linalg.norm(raw, axis=1, keepdims=True)
    norm[norm < 1e-8] = 1.0
    emb = raw / norm
    np.save(l2_path, emb)
    print(f"[INFO] 由 esm_emb.npy 生成 L2 归一化嵌入 -> {l2_path.name} "
          f"(范数 {np.linalg.norm(emb, axis=1).min():.3f}-{np.linalg.norm(emb, axis=1).max():.3f})")
    return emb


def main():
    emb = load_l2_embeddings()
    df = pd.read_csv(ANALYSIS / "umap_input.csv")
    assert len(emb) == len(df), f"长度不匹配 {len(emb)} vs {len(df)}"
    print(f"[INFO] 嵌入 {emb.shape}, 序列 {len(df)} 条")

    reducer = umap.UMAP(
        n_neighbors=30,
        min_dist=0.30,
        metric="cosine",
        random_state=42,
        n_components=2,
        verbose=True,
    )
    coords = reducer.fit_transform(emb)
    df["umap1"] = coords[:, 0]
    df["umap2"] = coords[:, 1]

    # 与长度的相关性 (验证长度是否主导)
    from scipy import stats
    r1 = stats.spearmanr(df["length"], df["umap1"])
    r2 = stats.spearmanr(df["length"], df["umap2"])
    print(f"\n[CHECK] 长度 vs UMAP1: ρ={r1.statistic:+.3f} (P={r1.pvalue:.1e})")
    print(f"[CHECK] 长度 vs UMAP2: ρ={r2.statistic:+.3f} (P={r2.pvalue:.1e})")

    # 新颖性: 我方候选到最近排行榜肽的距离
    lb = df[df["group"] == "leaderboard"][["umap1", "umap2"]].values
    ours = df[df["group"].str.startswith("ours")].copy()
    d = np.linalg.norm(ours[["umap1", "umap2"]].values[:, None, :] - lb[None, :, :], axis=2)
    df.loc[ours.index, "nn_dist_lb"] = d.min(axis=1)
    print(f"\n[CHECK] 我方候选到最近排行榜肽的 UMAP 距离: "
          f"中位 {np.median(d.min(axis=1)):.2f}, 范围 {d.min(axis=1).min():.2f}-{d.min(axis=1).max():.2f}")

    df.to_csv(ANALYSIS / "umap_coords.csv", index=False)
    print(f"\n[OK] 保存 {ANALYSIS}/umap_coords.csv")


if __name__ == "__main__":
    main()
