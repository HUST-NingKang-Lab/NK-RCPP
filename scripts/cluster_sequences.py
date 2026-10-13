#!/usr/bin/env python3
"""
纯 Python 序列聚类（替代 CD-HIT），按序列相似度去冗余
使用简单的贪婪聚类算法
"""
import argparse
from Bio import SeqIO
from collections import defaultdict


def kmer_similarity(seq1, seq2, k=3):
    """基于 k-mer 的 Jaccard 相似度"""
    if len(seq1) < k or len(seq2) < k:
        return 0.0
    kmers1 = set(str(seq1[i:i+k]) for i in range(len(seq1) - k + 1))
    kmers2 = set(str(seq2[i:i+k]) for i in range(len(seq2) - k + 1))
    intersection = len(kmers1 & kmers2)
    union = len(kmers1 | kmers2)
    return intersection / union if union > 0 else 0.0


def greedy_cluster(records, threshold):
    """贪婪聚类：按长度降序排序，每个cluster取第一个为代表"""
    # Sort by length descending (longer sequences as representatives)
    sorted_records = sorted(records, key=lambda r: len(r.seq), reverse=True)
    clusters = []
    representatives = []
    cluster_ids = {}

    for record in sorted_records:
        assigned = False
        for i, rep in enumerate(representatives):
            if kmer_similarity(record.seq, rep.seq) >= threshold:
                clusters[i].append(record)
                cluster_ids[record.id] = i
                assigned = True
                break
        if not assigned:
            clusters.append([record])
            representatives.append(record)
            cluster_ids[record.id] = len(representatives) - 1

    return clusters, representatives


def main():
    parser = argparse.ArgumentParser(description="Python-based sequence clustering")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output clustered FASTA file")
    parser.add_argument("--threshold", type=float, default=0.5,
                        help="Similarity threshold (default: 0.5, equivalent to -c 0.5)")
    parser.add_argument("--kmer", type=int, default=3,
                        help="K-mer size (default: 3)")
    args = parser.parse_args()

    print(f"[INFO] Reading sequences from {args.input}...")
    records = list(SeqIO.parse(args.input, "fasta"))
    print(f"[INFO] Loaded {len(records)} sequences")

    if len(records) == 0:
        print("[WARN] No sequences found")
        return

    # Use identity threshold = 1 - similarity_threshold
    # CD-HIT -c 0.5 means sequences with >= 50% identity are clustered
    # Our k-mer similarity: >= threshold means cluster together
    kmer_threshold = args.threshold

    print(f"[INFO] Clustering with k-mer similarity threshold {kmer_threshold}...")
    clusters, representatives = greedy_cluster(records, kmer_threshold)

    print(f"[INFO] {len(records)} sequences -> {len(clusters)} clusters")

    # Write representatives
    with open(args.output, "w") as f:
        for rep in representatives:
            SeqIO.write(rep, f, "fasta")
    print(f"[INFO] Wrote {len(representatives)} representatives to {args.output}")

    # Show cluster sizes
    sizes = [len(c) for c in clusters]
    print(f"[INFO] Cluster sizes: min={min(sizes)}, max={max(sizes)}, avg={sum(sizes)/len(sizes):.1f}")


if __name__ == "__main__":
    main()
