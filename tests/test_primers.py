import pytest
from ampliqc.core.primers import (
    KNOWN_PRIMERS,
    IUPAC_DICT,
    iupac_to_regex,
    reverse_complement,
    scan_primers,
)
from ampliqc.context.amplicon import AmpliconContextAnalyzer
from ampliqc.context.longreads import LongReadsContextAnalyzer


def test_iupac_to_regex():
    # 515F has Y (C or T), M (A or C) -> GTGYCAGCMGCCGCGGTAA
    pattern = iupac_to_regex("GTGYCAGCMGCCGCGGTAA")
    assert "[CT]" in pattern
    assert "[AC]" in pattern


def test_reverse_complement():
    seq = "AGAGTTTGATC"
    rc = reverse_complement(seq)
    assert rc == "GATCAAACTCT"


def test_scan_primers_16s_v4():
    # 515F sequence
    seqs = ["GTGCCAGCAGCCGCGGTAATACAGAGGATGCAAGCGTTATCCGG" + "A" * 100] * 50
    detected = scan_primers(seqs, min_pct=1.0)
    assert len(detected) > 0
    names = [p["primer_name"] for p in detected]
    assert any("515F" in n for n in names)


def test_scan_primers_ont_flanks():
    # ONT 16S flank sequence
    seqs = ["TTTCTGTTGGTGCTGATATTGCAGAGTTTGATCMTGGCTCAG" + "T" * 200] * 20
    detected = scan_primers(seqs, search_window=250, min_pct=1.0)
    names = [p["primer_name"] for p in detected]
    assert any("ONT-16S-Flank" in n for n in names)


def test_custom_primers_scanning():
    custom_primer = "ATCGATCGATCGATCG"
    seqs = [custom_primer + "ACGTACGTACGT"] * 30
    detected = scan_primers(seqs, custom_primers=[custom_primer], min_pct=1.0)
    names = [p["primer_name"] for p in detected]
    assert "Custom_Primer_1" in names


def test_longreads_primer_detection():
    analyzer = LongReadsContextAnalyzer(technology="ont")
    stats = {
        "total_reads": 100,
        "n50": 1500,
        "mean_read_length": 1200,
        "max_read_length": 3000,
        "min_read_length": 200,
        "read_qual_histogram": [0] * 20 + [100],
    }
    # ONT 16S Flanked sequence
    sample_seqs = ["TTTCTGTTGGTGCTGATATTGCAGAGTTTGATC" + "A" * 1000] * 100

    res = analyzer.evaluate(stats, sample_seqs)
    assert "detected_primers" in res["metadata"]
    assert len(res["metadata"]["detected_primers"]) > 0
    assert "ONT-16S-Flank" in res["metadata"]["detected_primers"][0]["primer_name"]

