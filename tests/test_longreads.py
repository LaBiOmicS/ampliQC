import pytest
from ampliqc.context.longreads import LongReadsContextAnalyzer


def test_longreads_evaluation():
    analyzer = LongReadsContextAnalyzer(technology="ont")

    stats = {
        "total_reads": 5000,
        "total_bases": 25000000,
        "mean_read_length": 5000.0,
        "min_read_length": 300,
        "max_read_length": 45000,
        "n50": 8500,
        "l50": 1200,
        "read_qual_histogram": [0]*10 + [500]*5 + [4000]*10 + [500]*36,  # >80% Q10+
    }

    sample_seqs = ["A" * 5000] * 100

    result = analyzer.evaluate(stats, sample_seqs)

    assert result["context"] == "Long Reads (ONT)"
    assert result["overall_status"] == "PASS"
    assert result["n50"] == 8500
