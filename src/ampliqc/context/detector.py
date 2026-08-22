"""
Sequencing context auto-detection engine.
"""

from typing import Dict, Any, List, Tuple
from ampliqc.core.primers import scan_primers
from ampliqc.context.amplicon import AmpliconContextAnalyzer
from ampliqc.context.base import ContextAnalyzer
from ampliqc.context.longreads import LongReadsContextAnalyzer


def detect_context(
    stats: Dict[str, Any], sample_sequences: List[str]
) -> Tuple[str, ContextAnalyzer]:
    """
    Automatically detects Amplicon sequencing context (Short-Read vs Long-Read 16S/18S/ITS/Functional).
    """
    mean_len = stats.get("mean_read_length", 0)
    max_len = stats.get("max_read_length", 0)
    is_long_read = mean_len > 800 or max_len > 2500

    if not sample_sequences:
        target_default = "amplicon_longread" if is_long_read else "amplicon_shortread"
        return target_default, AmpliconContextAnalyzer(target_region=target_default, is_long_read=is_long_read)

    # 1. Scan for known amplicon primers and technology flanks
    detected = scan_primers(
        sequences=sample_sequences,
        search_window=250 if is_long_read else 120,
        min_pct=3.0,
    )
    primer_hits = {p["primer_name"]: p["match_pct"] for p in detected}

    # 2. If primers detected -> Amplicon Strategy
    if primer_hits:
        top_primer = max(primer_hits, key=primer_hits.get)
        if "16S" in top_primer:
            ctx_name = "amplicon_16s"
        elif "ITS" in top_primer:
            ctx_name = "amplicon_its"
        elif "18S" in top_primer:
            ctx_name = "amplicon_18s"
        else:
            ctx_name = "amplicon_generic"
            
        return ctx_name, AmpliconContextAnalyzer(target_region=ctx_name, is_long_read=is_long_read)

    # 3. Default Amplicon fallback for short and long reads
    if is_long_read:
        return "amplicon_longread", AmpliconContextAnalyzer(target_region="full_length_amplicon", is_long_read=True)
    else:
        return "amplicon_shortread", AmpliconContextAnalyzer(target_region="shortread_amplicon", is_long_read=False)

