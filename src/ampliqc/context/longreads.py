"""
Long Reads Context Analyzer for Oxford Nanopore (ONT) and PacBio HiFi / MAS-Seq.
"""

from typing import Dict, Any, List
from ampliqc.context.base import ContextAnalyzer, QualityCheckResult, Status
from ampliqc.core.primers import scan_primers, PRIMER_REGION_MAP


class LongReadsContextAnalyzer(ContextAnalyzer):
    """
    Quality Control evaluator specialized for Long Reads (Oxford Nanopore ONT / PacBio HiFi).
    """

    def __init__(self, technology: str = "long_reads"):
        super().__init__(context_name=f"Long Reads ({technology.upper()})")
        self.technology = technology

    def evaluate(self, stats: Dict[str, Any], sample_sequences: List[str]) -> Dict[str, Any]:
        results = []

        # 1. Read Length Distribution & N50 Evaluation
        n50_res = self._check_length_and_n50(stats)
        results.append(n50_res)

        # 2. Quality Score Profile (Adapted for Long Reads Q10-Q20+)
        qual_res = self._check_long_read_quality(stats)
        results.append(qual_res)

        # 3. Short Fragment Contamination Check (< 500bp)
        short_frag_res = self._check_short_fragments(stats)
        results.append(short_frag_res)

        # 4. Long Read Primer & Adapter Scan (ONT Flanks, PacBio Adapters, Full-Length 16S/18S/ITS)
        primer_res = self._check_longread_primers(sample_sequences)
        results.append(primer_res)

        # Overall Status
        statuses = [r.status for r in results]
        if Status.FAIL in statuses:
            overall = Status.FAIL
        elif Status.WARN in statuses:
            overall = Status.WARN
        else:
            overall = Status.PASS

        detected_primers = primer_res.metrics.get("detected_primers", [])
        primary_primer_str = "None detected (Clean or already demultiplexed)"
        has_untrimmed = False

        if detected_primers:
            top_primer = max(detected_primers, key=lambda x: x["match_pct"])
            primary_primer_str = f"{top_primer['primer_name']} [{top_primer['sequence']}] ({top_primer['match_pct']}% match)"
            if top_primer["match_pct"] > 5.0:
                has_untrimmed = True

        trimming_status = "UNTRIMMED (Chopper / pycoQC / Cutadapt Trimming Needed)" if has_untrimmed else "TRIMMED / CLEAN"

        metadata = {
            "sequencing_technology": f"Long-Read ({self.technology.upper()})",
            "context_type": "Long-Read Sequencing",
            "target_region": "Full Length Genomes / Structural Variants / Full-Length Amplicons",
            "detected_primers": detected_primers,
            "primary_primer": primary_primer_str,
            "trimming_status": trimming_status,
            "recommended_pipeline": "Dorado / Chopper ➔ Flye / Minimap2 / pycoQC",
        }

        return {
            "context": self.context_name,
            "overall_status": overall.value,
            "checks": [r.to_dict() for r in results],
            "n50": stats.get("n50", 0),
            "l50": stats.get("l50", 0),
            "primer_metrics": primer_res.metrics,
            "metadata": metadata,
        }

    def _check_longread_primers(self, sequences: List[str]) -> QualityCheckResult:
        """
        Scans long reads (ONT / PacBio) for sequencing adapters, flanking barcodes, and PCR primers.
        """
        detected_primers = scan_primers(
            sequences=sequences,
            search_window=250,
            min_pct=1.0
        )

        if not detected_primers:
            status = Status.PASS
            msg = "No residual long-read primers or flanking adapters detected."
            reasoning = "Reads appear cleanly trimmed or demultiplexed."
        else:
            top_primer = max(detected_primers, key=lambda x: x["match_pct"])
            if top_primer["match_pct"] > 25.0:
                status = Status.WARN
                msg = f"High residual long-read adapter/primer content ({top_primer['primer_name']} in {top_primer['match_pct']}% of reads)."
                reasoning = "High residual adapter content in long reads may impair de novo assembly or alignment alignment quality. Recommend Chopper / Cutadapt."
            else:
                status = Status.PASS
                msg = f"Minor long-read primer/adapter residue ({top_primer['primer_name']} in {top_primer['match_pct']}% of reads)."
                reasoning = "Low residual adapter rate is within normal tolerance."

        return QualityCheckResult(
            name="Long-Read Primer & Flank Adapter Scan",
            status=status,
            message=msg,
            metrics={"detected_primers": detected_primers},
            context_reasoning=reasoning
        )

    def _check_length_and_n50(self, stats: Dict[str, Any]) -> QualityCheckResult:
        n50 = stats.get("n50", 0)
        mean_len = stats.get("mean_read_length", 0)
        max_len = stats.get("max_read_length", 0)
        min_len = stats.get("min_read_length", 0)

        if n50 >= 1000:
            status = Status.PASS
            msg = f"N50 is {n50:,} bp (Mean: {mean_len:,.1f} bp, Max: {max_len:,} bp)."
            reasoning = "N50 read length profile is excellent for long-read assembly, full-length 16S/ITS, or structural variant analysis."
        else:
            status = Status.WARN
            msg = f"Low N50 read length ({n50:,} bp)."
            reasoning = "N50 is below 1,000 bp. Check DNA extraction/fragmentation protocol or size selection."

        return QualityCheckResult(
            name="Read Length N50 & Span",
            status=status,
            message=msg,
            metrics={"n50": n50, "mean": mean_len, "max": max_len, "min": min_len},
            context_reasoning=reasoning
        )

    def _check_long_read_quality(self, stats: Dict[str, Any]) -> QualityCheckResult:
        read_qual_hist = stats.get("read_qual_histogram", [])
        if not read_qual_hist:
            return QualityCheckResult(
                name="Long Read Quality Profile",
                status=Status.WARN,
                message="No read quality histogram data.",
                metrics={},
                context_reasoning="Missing quality data."
            )

        total_reads = max(1, stats.get("total_reads", 1))
        # Reads with Phred >= Q10 (90% accuracy) and Q15 (97% accuracy)
        q10_reads = sum(read_qual_hist[10:])
        q15_reads = sum(read_qual_hist[15:])

        q10_pct = round((q10_reads / total_reads) * 100, 1)
        q15_pct = round((q15_reads / total_reads) * 100, 1)

        if q10_pct >= 80.0:
            status = Status.PASS
            msg = f"{q10_pct}% of reads are ≥ Q10 ({q15_pct}% ≥ Q15)."
            reasoning = (
                "Long read accuracy profile is high. FastQC often marks Q10-Q20 average scores as low, but Q10-Q15+ "
                "is standard and high quality for Oxford Nanopore R10.4.1 and PacBio HiFi."
            )
        elif q10_pct >= 50.0:
            status = Status.WARN
            msg = f"{q10_pct}% of reads are ≥ Q10."
            reasoning = "Moderate read quality. Consider filtering low Q-score reads (e.g. Chopper / NanoFilt q < 10)."
        else:
            status = Status.FAIL
            msg = f"Low read quality: Only {q10_pct}% of reads are ≥ Q10."
            reasoning = "High proportion of sub-Q10 reads. Pore degradation or flowcell issues suspected."

        return QualityCheckResult(
            name="Long Read Accuracy (Q10/Q15 Profile)",
            status=status,
            message=msg,
            metrics={"q10_percentage": q10_pct, "q15_percentage": q15_pct},
            context_reasoning=reasoning
        )

    def _check_short_fragments(self, stats: Dict[str, Any]) -> QualityCheckResult:
        min_len = stats.get("min_read_length", 0)
        mean_len = stats.get("mean_read_length", 0)

        if min_len < 200 and mean_len > 2000:
            status = Status.WARN
            msg = f"Short fragment contamination detected (Min length: {min_len} bp)."
            reasoning = (
                "Presence of very short reads (< 200 bp) in long-read sequencing consumes flowcell pores inefficiently. "
                "Recommend applying size selection or filtering short reads before assembly."
            )
        else:
            status = Status.PASS
            msg = f"Fragment length distribution is consistent (Min: {min_len} bp)."
            reasoning = "No significant pore-clogging short fragment contamination detected."

        return QualityCheckResult(
            name="Short Fragment Contamination",
            status=status,
            message=msg,
            metrics={"min_length": min_len},
            context_reasoning=reasoning
        )
