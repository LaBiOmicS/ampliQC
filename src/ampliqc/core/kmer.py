"""
K-mer analysis module for richness, overrepresented k-mers, and bias detection.
"""

from collections import Counter
from typing import Dict, List, Any


def analyze_kmers(
    sequences: List[str], k: int = 7, top_n: int = 20
) -> Dict[str, Any]:
    """
    Computes k-mer frequency, richness, and top overrepresented k-mers.
    """
    kmer_counts = Counter()
    total_kmers = 0

    for seq in sequences:
        seq_len = len(seq)
        if seq_len < k:
            continue
        for i in range(seq_len - k + 1):
            kmer = seq[i : i + k]
            if "N" not in kmer:
                kmer_counts[kmer] += 1
                total_kmers += 1

    unique_kmers = len(kmer_counts)
    max_possible = 4**k
    richness_ratio = round(unique_kmers / max_possible, 4) if max_possible > 0 else 0

    top_kmers = []
    for kmer, count in kmer_counts.most_common(top_n):
        pct = round(100.0 * count / max(1, total_kmers), 3)
        top_kmers.append({"kmer": kmer, "count": count, "percentage": pct})

    return {
        "k_size": k,
        "total_kmers": total_kmers,
        "unique_kmers": unique_kmers,
        "kmer_richness_ratio": richness_ratio,
        "top_kmers": top_kmers,
    }
