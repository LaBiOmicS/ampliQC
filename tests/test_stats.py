import pytest
from ampliqc.core.parser import FastqRecord
from ampliqc.core.stats import FastqStatsAccumulator


def test_stats_accumulator():
    acc = FastqStatsAccumulator()
    rec1 = FastqRecord("@seq1", "CCTACGGGAGGCAGCAGTAGGGAATATTGG", "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII")
    rec2 = FastqRecord("@seq2", "CCTACGGGAGGCAGCAGTAGGGAATATTGG", "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII")

    acc.process_record(rec1)
    acc.process_record(rec2)

    summary = acc.get_summary()

    assert summary["total_reads"] == 2
    assert summary["total_bases"] == 60
    assert summary["mean_read_length"] == 30.0
    assert summary["duplication_rate"] == 50.0
    assert summary["overall_gc_content"] > 0
