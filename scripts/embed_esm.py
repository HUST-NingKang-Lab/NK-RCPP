#!/usr/bin/env python3
"""
ESM2-650M 平均池化编码
模型与 RAPiDock 完全一致: esm2_t33_650M_UR50D, 第33层 (repr_layers=[33])
池化: mean pooling over residue tokens (排除 BOS/EOS)

需在有 GPU 的机器上运行 (conda activate RAPiDock)
用法: python embed_esm.py --input analysis/umap_input.csv --output analysis/esm_emb.npy
"""
import argparse
import numpy as np
import pandas as pd
import torch
import esm
from pathlib import Path


def mean_pool(rep, n_res):
    """rep: (1, L+2, D) 含 BOS/EOS -> 对残基位置取均值"""
    return rep[0, 1:n_res + 1].mean(dim=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--model", default="esm2_t33_650M_UR50D")
    ap.add_argument("--batch_tokens", type=int, default=4096)
    args = ap.parse_args()

    df = pd.read_csv(args.input)
    seqs = df["seq"].astype(str).tolist()
    print(f"[INFO] {len(seqs)} 条序列, 长度 {min(map(len,seqs))}-{max(map(len,seqs))}")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[INFO] device = {device}")

    model, alphabet = esm.pretrained.load_model_and_alphabet(args.model)
    model.eval().to(device)
    print(f"[INFO] 模型加载完成: {args.model}, 层数 {model.num_layers}, 维度 {model.embed_dim}")

    batch_converter = alphabet.get_batch_converter()
    repr_layer = model.num_layers  # 与 RAPiDock 的 repr_layers=[33] 一致

    embs = np.zeros((len(seqs), model.embed_dim), dtype=np.float32)
    # 按 token 预算分批
    order = np.argsort([len(s) for s in seqs])  # 长度排序, 减少 padding
    batches, cur, cur_tok = [], [], 0
    for idx in order:
        t = len(seqs[idx]) + 2
        if cur and cur_tok + t > args.batch_tokens:
            batches.append(cur); cur, cur_tok = [], 0
        cur.append(idx); cur_tok += t
    if cur:
        batches.append(cur)

    with torch.no_grad():
        for bi, batch_idx in enumerate(batches):
            data = [(str(i), seqs[i]) for i in batch_idx]
            _, _, toks = batch_converter(data)
            toks = toks.to(device)
            out = model(toks, repr_layers=[repr_layer], return_contacts=False)
            reps = out["representations"][repr_layer]
            for j, i in enumerate(batch_idx):
                embs[i] = mean_pool(reps[j], len(seqs[i])).cpu().numpy()
            print(f"  batch {bi+1}/{len(batches)} ({len(batch_idx)} seqs)")

    # sanity checks
    norms = np.linalg.norm(embs, axis=1)
    print(f"[CHECK] embedding 范数: min={norms.min():.2f} max={norms.max():.2f} "
          f"nan={np.isnan(embs).sum()}")
    # 同一序列应完全一致 (用重复序列验证)
    dup = df[df.duplicated("seq", keep=False)]
    if len(dup):
        s = dup["seq"].iloc[0]
        ii = df.index[df["seq"] == s].tolist()
        d = np.linalg.norm(embs[ii[0]] - embs[ii[1]])
        print(f"[CHECK] 重复序列({s}) 编码距离 = {d:.6f} (应≈0)")

    np.save(args.output, embs)
    print(f"[OK] 保存 {args.output}  shape={embs.shape}")


if __name__ == "__main__":
    main()
