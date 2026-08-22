"""
Amplicon-specific context analyzer (16S, 18S, ITS).
"""

import re
from typing import Dict, Any, List, Optional, Tuple
from ampliqc.context.base import ContextAnalyzer, QualityCheckResult, Status

from ampliqc.reports.helpers import calculate_dada2_recommendation, generate_cutadapt_command
from ampliqc.core.primers import (
    IUPAC_DICT,
    KNOWN_PRIMERS,
    PRIMER_REGION_MAP,
    iupac_to_regex,
    reverse_complement,
    scan_primers,
)

# Target region length ranges (bp)
EXPECTED_LENGTHS = {
    "16s_v3v4": (400, 470),
    "16s_v4": (240, 290),
    "its1": (180, 400),
    "its2": (200, 450),
    "18s_v9": (120, 180),
    "generic_amplicon": (100, 600),
}


class AmpliconContextAnalyzer(ContextAnalyzer):
    """
    Quality Control evaluator specialized for 16S, 18S, ITS, and marker gene amplicon sequencing (Short and Long Reads).
    """

    def __init__(
        self,
        target_region: str = "auto",
        custom_primers: Optional[List[str]] = None,
        is_long_read: bool = False,
    ):
        super().__init__(context_name=f"Amplicon ({target_region})")
        self.target_region = target_region
        self.custom_primers = custom_primers or []
        self.is_long_read = is_long_read

    def evaluate(self, stats: Dict[str, Any], sample_sequences: List[str]) -> Dict[str, Any]:
        results = []

        # 1. Primer Match & Efficiency Analysis
        primer_res = self._check_primers(sample_sequences)
        results.append(primer_res)

        # 2. Sequence Duplication Re-evaluation (Context-aware)
        dup_res = self._check_duplication(stats)
        results.append(dup_res)

        # 3. Per-base Content Re-evaluation (Primer-aware bias)
        content_res = self._check_per_base_content(stats, primer_res.metrics)
        results.append(content_res)

        # 4. Read Length vs Expected Amplicon Target Length
        len_res = self._check_read_lengths(stats)
        results.append(len_res)

        # Overall Status Determination
        statuses = [r.status for r in results]
        if Status.FAIL in statuses:
            overall = Status.FAIL
        elif Status.WARN in statuses:
            overall = Status.WARN
        else:
            overall = Status.PASS

        detected_primers = primer_res.metrics.get("detected_primers", [])
        
        # Infer target region from detected primers if auto or general
        inferred_region = self.target_region
        primary_primer_str = "None detected (Already trimmed or clean)"

        if detected_primers:
            fwd_primers = [p for p in detected_primers if "F" in p["primer_name"] or "Fwd" in p["primer_name"]]
            rev_primers = [p for p in detected_primers if "R" in p["primer_name"] or "Rev" in p["primer_name"]]

            top_fwd = max(fwd_primers, key=lambda x: x["match_pct"]) if fwd_primers else None
            top_rev = max(rev_primers, key=lambda x: x["match_pct"]) if rev_primers else None

            if top_fwd and top_rev:
                primary_primer_str = f"Fwd: {top_fwd['primer_name']} ({top_fwd['match_pct']}%) | Rev: {top_rev['primer_name']} ({top_rev['match_pct']}%)"
                primary = top_fwd
            elif top_fwd:
                primary_primer_str = f"Fwd: {top_fwd['primer_name']} [{top_fwd['sequence']}] ({top_fwd['match_pct']}% match)"
                primary = top_fwd
            else:
                top_primer = max(detected_primers, key=lambda x: x["match_pct"])
                primary_primer_str = f"{top_primer['primer_name']} [{top_primer['sequence']}] ({top_primer['match_pct']}% match)"
                primary = top_primer

            for p_code, r_name in PRIMER_REGION_MAP.items():
                if p_code in primary["primer_name"]:
                    inferred_region = r_name
                    break

        mean_len = stats.get("mean_read_length", 0)
        is_long = self.is_long_read or mean_len > 800

        if inferred_region in ["auto", "amplicon_16s"]:
            inferred_region = "16S rRNA (Full-Length ~1.5kb)" if is_long else "16S rRNA"
        elif inferred_region == "amplicon_its":
            inferred_region = "ITS (Full-Length)" if is_long else "ITS (Fungi)"
        elif inferred_region == "amplicon_18s":
            inferred_region = "18S rRNA (Full-Length)" if is_long else "18S rRNA (Eukaryotes)"

        has_untrimmed = any(p["match_pct"] > 5.0 for p in detected_primers)
        trimming_status = "UNTRIMMED (Trimming Required)" if has_untrimmed else "TRIMMED / CLEAN"

        if is_long:
            seq_tech = f"Long-Read (PacBio HiFi / ONT, ~{mean_len:.0f} bp)"
            ctx_label = "Amplicon Sequencing (Long-Read)"
            rec_pipe = "Cutadapt / Chopper ➔ DADA2 (Long-Read) / PB-16S-nf / QIIME 2 (Full-Length ASVs)"
        else:
            seq_tech = f"Illumina Short-Read (~{mean_len:.0f} bp)"
            ctx_label = "Amplicon Sequencing (Short-Read)"
            rec_pipe = "Cutadapt (Primer Trimming) ➔ DADA2 / Deblur (ASV Clustering) ➔ QIIME 2 / SILVA / UNITE"

        dada2_rec = calculate_dada2_recommendation(stats, inferred_region, is_paired=not is_long)
        cutadapt_cmd = generate_cutadapt_command(
            sample_name="sample",
            r1_path="reads_R1.fastq.gz",
            r2_path="reads_R2.fastq.gz" if not is_long else None,
            detected_primers=detected_primers,
            is_long_read=is_long,
        )

        metadata = {
            "sequencing_technology": seq_tech,
            "context_type": ctx_label,
            "target_region": inferred_region,
            "detected_primers": detected_primers,
            "primary_primer": primary_primer_str,
            "trimming_status": trimming_status,
            "recommended_pipeline": rec_pipe,
            "dada2_recommendation": dada2_rec["recommendation_text"],
            "dada2_code": dada2_rec["dada2_code"],
            "cutadapt_command": cutadapt_cmd,
        }

        return {
            "context": self.context_name,
            "overall_status": overall.value,
            "checks": [r.to_dict() for r in results],
            "primer_metrics": primer_res.metrics,
            "metadata": metadata,
        }

    def _check_primers(self, sequences: List[str]) -> QualityCheckResult:
        """
        Scans reads for known amplicon primers and evaluates trimming status and orientation.
        """
        if not sequences:
            return QualityCheckResult(
                name="Primer Presence & Trimming",
                status=Status.WARN,
                message="No sequences available to evaluate primer presence.",
                metrics={},
                context_reasoning="Empty sequence sample."
            )

        detected_primers = scan_primers(
            sequences=sequences,
            custom_primers=self.custom_primers,
            search_window=200 if self.is_long_read else 150,
            min_pct=1.0,
        )

        if not detected_primers:
            status = Status.PASS
            msg = "No untrimmed primers detected in sequence heads/tails (Primers already trimmed or clean)."
            reasoning = "Amplicon reads appear to be cleanly trimmed of PCR primers."
        else:
            top_primer = max(detected_primers, key=lambda x: x["match_pct"])
            if top_primer["match_pct"] > 30.0:
                status = Status.WARN
                msg = f"High proportion of untrimmed primers detected ({top_primer['primer_name']} in {top_primer['match_pct']}% of reads)."
                reasoning = (
                    "Untrimmed primers remaining in amplicon reads can distort downstream OTU/ASV feature table generation (e.g. DADA2 / Deblur). "
                    "Recommend running primer trimming (e.g. Cutadapt) before ASV clustering."
                )
            else:
                status = Status.PASS
                msg = f"Minor primer residue detected ({top_primer['primer_name']} in {top_primer['match_pct']}% of reads)."
                reasoning = "Low residual primer presence is acceptable prior to trimming."

        return QualityCheckResult(
            name="Primer Detection & Trimming Status",
            status=status,
            message=msg,
            metrics={"detected_primers": detected_primers},
            context_reasoning=reasoning
        )

    def _check_duplication(self, stats: Dict[str, Any]) -> QualityCheckResult:
        """
        Evaluates sequence duplication in an amplicon-aware context.
        Standard FastQC flags > 50% duplication as FAIL.
        In amplicon sequencing, high duplication is expected and desired because reads originate from specific amplified loci.
        """
        dup_rate = stats.get("duplication_rate", 0.0)
        total_reads = stats.get("total_reads", 0)

        # In amplicons, dup_rate > 80% is common and indicates dominant ASVs/OTUs
        if dup_rate > 98.0:
            status = Status.WARN
            msg = f"Extremely high duplication rate ({dup_rate}%)."
            reasoning = (
                "While high duplication is natural in amplicon libraries due to locus targeting, >98% duplication may indicate low library complexity, "
                "monoclonal PCR explosion, or index hopping."
            )
        else:
            status = Status.PASS
            msg = f"Duplication rate is {dup_rate}% (Normal/Expected for amplicon locus amplification)."
            reasoning = (
                "Standard FastQC flags high duplication rates as a failure. However, amplicon sequencing (16S/18S/ITS) inherently targets specific genomic regions, "
                "making high sequence duplication an expected biological characteristic rather than a library error."
            )

        return QualityCheckResult(
            name="Sequence Duplication (Amplicon Context)",
            status=status,
            message=msg,
            metrics={"duplication_rate": dup_rate, "total_reads": total_reads},
            context_reasoning=reasoning
        )

    def _check_per_base_content(self, stats: Dict[str, Any], primer_metrics: Dict[str, Any]) -> QualityCheckResult:
        """
        Evaluates per-base sequence content considering conserved amplicon primer/region signatures.
        """
        content = stats.get("per_base_content", {})
        if not content or not content.get("A"):
            return QualityCheckResult(
                name="Per-base Sequence Content",
                status=Status.PASS,
                message="No base content metrics.",
                metrics={},
                context_reasoning="Insufficient data."
            )

        # Check variance in positions 25 to end (downstream of primers)
        seq_len = len(content["A"])
        body_start = min(30, seq_len // 4)
        
        a_body = content["A"][body_start:]
        c_body = content["C"][body_start:]
        g_body = content["G"][body_start:]
        t_body = content["T"][body_start:]

        # Calculate difference between A/T and C/G in read body
        if a_body and c_body:
            max_dev = max(
                max(abs(a - t) for a, t in zip(a_body, t_body)),
                max(abs(c - g) for c, g in zip(c_body, g_body))
            )
        else:
            max_dev = 0

        has_untrimmed_primers = len(primer_metrics.get("detected_primers", [])) > 0

        if max_dev > 35.0:
            status = Status.WARN
            msg = "Significant base content imbalance in read body."
            reasoning = "Large base composition divergence outside primer region."
        else:
            status = Status.PASS
            msg = "Per-base sequence content is within expected amplicon composition limits."
            reasoning = (
                "FastQC often flags position 1-25 nucleotide bias as a failure in amplicons due to conserved primer binding sites. "
                "ampliQC accounts for primer composition and verifies stability across the read body."
            )

        return QualityCheckResult(
            name="Per-base Sequence Content (Amplicon Context)",
            status=status,
            message=msg,
            metrics={"max_body_deviation": round(max_dev, 2)},
            context_reasoning=reasoning
        )

    def _check_read_lengths(self, stats: Dict[str, Any]) -> QualityCheckResult:
        mean_len = stats.get("mean_read_length", 0)
        min_len = stats.get("min_read_length", 0)
        max_len = stats.get("max_read_length", 0)

        status = Status.PASS
        msg = f"Read length distribution: Mean {mean_len}bp (Min: {min_len}bp, Max: {max_len}bp)."
        reasoning = "Read lengths are suitable for amplicon processing."

        return QualityCheckResult(
            name="Read Length Profile",
            status=status,
            message=msg,
            metrics={"mean": mean_len, "min": min_len, "max": max_len},
            context_reasoning=reasoning
        )

