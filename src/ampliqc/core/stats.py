"""
Core QC statistics accumulator using fast numpy vectorized math.
"""

from collections import Counter
import numpy as np
from typing import Dict, Any, List, Optional
from ampliqc.core.parser import read_fastq, FastqRecord
from ampliqc.core.adapters import scan_adapter_content


class FastqStatsAccumulator:
    """
    Accumulates metrics across streamed FASTQ records.
    """

    def __init__(self, sample_limit: Optional[int] = None):
        self.sample_limit = sample_limit
        self.total_reads = 0
        self.total_bases = 0
        self.read_lengths = []
        
        # Max expected read length dynamic array sizing
        self._max_seen_len = 0
        
        # Per-base ASCII quality sums and counts: matrix [max_len, 128]
        self._qual_matrix: Optional[np.ndarray] = None
        
        # Per-base base counts: matrix [max_len, 5] (A, C, G, T, N)
        self._base_matrix: Optional[np.ndarray] = None
        
        # Per-read average quality distribution (histogram 0 to 60)
        self.read_qual_hist = np.zeros(61, dtype=np.int64)
        
        # Per-read GC content distribution (histogram 0 to 100%)
        self.gc_hist = np.zeros(101, dtype=np.int64)
        
        # Exact sequence counter for duplicate detection
        self.sequence_counts = Counter()
        self.sampled_sequences = 0

    def process_record(self, record: FastqRecord):
        self.total_reads += 1
        seq = record.sequence.upper()
        qual_str = record.qualities
        seq_len = len(seq)
        
        self.total_bases += seq_len
        self.read_lengths.append(seq_len)
        
        if seq_len > self._max_seen_len:
            self._resize_matrices(seq_len)
            self._max_seen_len = seq_len

        # Convert ASCII qualities to Phred scores (Phred+33)
        quals = np.frombuffer(qual_str.encode('latin1'), dtype=np.uint8) - 33
        
        # Vectorized per-base quality accumulation
        for pos, q in enumerate(quals):
            self._qual_matrix[pos, min(q, 60)] += 1
            
        mean_q = int(round(np.mean(quals))) if len(quals) > 0 else 0
        mean_q = max(0, min(60, mean_q))
        self.read_qual_hist[mean_q] += 1
        
        # Base counts & GC content
        gc_count = 0
        for pos, char in enumerate(seq):
            if char == 'A':
                self._base_matrix[pos, 0] += 1
            elif char == 'C':
                self._base_matrix[pos, 1] += 1
                gc_count += 1
            elif char == 'G':
                self._base_matrix[pos, 2] += 1
                gc_count += 1
            elif char == 'T':
                self._base_matrix[pos, 3] += 1
            else:  # N or non-ACGT
                self._base_matrix[pos, 4] += 1
                
        gc_pct = int(round((gc_count / seq_len) * 100)) if seq_len > 0 else 0
        self.gc_hist[min(100, gc_pct)] += 1

        # Track sequences for duplication analysis (sample first 200,000 to save memory)
        if self.sampled_sequences < 200000:
            self.sequence_counts[seq] += 1
            self.sampled_sequences += 1

    def _resize_matrices(self, new_max_len: int):
        new_max_len = max(new_max_len, self._max_seen_len + 50)
        
        if self._qual_matrix is None:
            self._qual_matrix = np.zeros((new_max_len, 61), dtype=np.int64)
            self._base_matrix = np.zeros((new_max_len, 5), dtype=np.int64)
        else:
            old_len = self._qual_matrix.shape[0]
            if new_max_len > old_len:
                new_qual = np.zeros((new_max_len, 61), dtype=np.int64)
                new_base = np.zeros((new_max_len, 5), dtype=np.int64)
                new_qual[:old_len, :] = self._qual_matrix
                new_base[:old_len, :] = self._base_matrix
                self._qual_matrix = new_qual
                self._base_matrix = new_base

    def get_summary(self) -> Dict[str, Any]:
        """
        Calculates and returns structured QC summary dictionary.
        """
        if self.total_reads == 0:
            return {"total_reads": 0, "error": "Empty FASTQ file"}
            
        max_pos = self._max_seen_len
        qual_mat = self._qual_matrix[:max_pos, :]
        base_mat = self._base_matrix[:max_pos, :]
        
        # Calculate per-base quality quantiles & means
        per_base_qual_mean = []
        per_base_qual_median = []
        per_base_qual_q25 = []
        per_base_qual_q75 = []
        
        for pos in range(max_pos):
            counts = qual_mat[pos, :]
            total_at_pos = np.sum(counts)
            if total_at_pos == 0:
                per_base_qual_mean.append(0.0)
                per_base_qual_median.append(0)
                per_base_qual_q25.append(0)
                per_base_qual_q75.append(0)
                continue
            
            scores = np.arange(61)
            mean_score = float(np.sum(scores * counts) / total_at_pos)
            per_base_qual_mean.append(round(mean_score, 2))
            
            # Cumulative distribution for percentiles
            cdf = np.cumsum(counts) / total_at_pos
            med = int(np.searchsorted(cdf, 0.50))
            q25 = int(np.searchsorted(cdf, 0.25))
            q75 = int(np.searchsorted(cdf, 0.75))
            
            per_base_qual_median.append(med)
            per_base_qual_q25.append(q25)
            per_base_qual_q75.append(q75)
            
        # Calculate per-base base content percentages
        per_base_content = {"A": [], "C": [], "G": [], "T": [], "N": []}
        for pos in range(max_pos):
            b_counts = base_mat[pos, :]
            total_b = np.sum(b_counts)
            if total_b == 0:
                for k in per_base_content:
                    per_base_content[k].append(0.0)
            else:
                per_base_content["A"].append(round(float(b_counts[0] / total_b * 100), 2))
                per_base_content["C"].append(round(float(b_counts[1] / total_b * 100), 2))
                per_base_content["G"].append(round(float(b_counts[2] / total_b * 100), 2))
                per_base_content["T"].append(round(float(b_counts[3] / total_b * 100), 2))
                per_base_content["N"].append(round(float(b_counts[4] / total_b * 100), 2))

        # Overall GC %
        total_acgt = np.sum(base_mat[:, :4])
        gc_bases = np.sum(base_mat[:, 1]) + np.sum(base_mat[:, 2])
        overall_gc = round(float((gc_bases / total_acgt * 100)), 2) if total_acgt > 0 else 0.0

        # Duplication & overrepresented sequence calculation
        unique_seqs = len(self.sequence_counts)
        dup_rate = round(100.0 * (1.0 - (unique_seqs / max(1, self.sampled_sequences))), 2)
        
        overrepresented = []
        for seq, count in self.sequence_counts.most_common(20):
            pct = round(100.0 * count / max(1, self.sampled_sequences), 3)
            if pct >= 0.1:  # Sequences representing >= 0.1% of reads
                overrepresented.append({
                    "sequence": seq,
                    "count": count,
                    "percentage": pct,
                    "length": len(seq)
                })

        # Calculate adapter profiles across positions
        sample_seq_list = list(self.sequence_counts.keys())[:50000]
        adapter_profiles = scan_adapter_content(sample_seq_list, max_pos)

        # Calculate N50 and L50 (essential for long reads ONT/PacBio)
        n50 = 0
        l50 = 0
        median_len = 0
        if self.read_lengths:
            sorted_lens = np.sort(self.read_lengths)[::-1]
            cum_sum = np.cumsum(sorted_lens)
            target = self.total_bases / 2.0
            idx = int(np.searchsorted(cum_sum, target))
            if idx < len(sorted_lens):
                n50 = int(sorted_lens[idx])
                l50 = idx + 1
            median_len = int(np.median(self.read_lengths))

        return {
            "total_reads": self.total_reads,
            "total_bases": self.total_bases,
            "mean_read_length": round(float(np.mean(self.read_lengths)), 1) if self.read_lengths else 0,
            "median_read_length": median_len,
            "min_read_length": int(np.min(self.read_lengths)) if self.read_lengths else 0,
            "max_read_length": int(np.max(self.read_lengths)) if self.read_lengths else 0,
            "n50": n50,
            "l50": l50,
            "overall_gc_content": overall_gc,
            "duplication_rate": dup_rate,
            "per_base_qual_mean": per_base_qual_mean,
            "per_base_qual_median": per_base_qual_median,
            "per_base_qual_q25": per_base_qual_q25,
            "per_base_qual_q75": per_base_qual_q75,
            "read_qual_histogram": self.read_qual_hist.tolist(),
            "gc_content_histogram": self.gc_hist.tolist(),
            "per_base_content": per_base_content,
            "adapter_content": adapter_profiles,
            "overrepresented_sequences": overrepresented,
            "sequence_counts": self.sequence_counts,
        }
