"""
Standard Illumina, Nextera, Small RNA, and Poly-A adapter scanner.
"""

from typing import Dict, List, Any

# Standard sequencing adapters
KNOWN_ADAPTERS = {
    "Illumina Universal Adapter": "AGATCGGAAGAG",
    "Illumina Small RNA 3' Adapter": "TGGGAATCTCGGG",
    "Illumina Small RNA 5' Adapter": "GUUCAGAGUUCUACAGUCCGACGAUC",
    "Nextera Transposase Sequence": "CTGTCTCTTATACACATCT",
    "TruSeq Read 1 Adapter": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    "TruSeq Read 2 Adapter": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    "Poly-A Tail": "AAAAAAAAAAAAAAAAAAAA",
}


def scan_adapter_content(
    sequences: List[str], max_len: int
) -> Dict[str, List[float]]:
    """
    Computes per-position adapter cumulative percentage across reads.
    """
    n_reads = max(1, len(sequences))
    adapter_profiles = {}

    for name, adapter_seq in KNOWN_ADAPTERS.items():
        k = len(adapter_seq)
        pos_counts = [0] * max_len

        for seq in sequences:
            seq_len = len(seq)
            # Find adapter occurrence in sequence
            idx = seq.find(adapter_seq[:8])  # Match seed 8-mer
            if idx != -1:
                for p in range(idx, seq_len):
                    pos_counts[p] += 1

        pcts = [round(100.0 * c / n_reads, 2) for c in pos_counts]
        adapter_profiles[name] = pcts

    return adapter_profiles
