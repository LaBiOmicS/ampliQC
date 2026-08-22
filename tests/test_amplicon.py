import pytest
from ampliqc.context.amplicon import AmpliconContextAnalyzer


def test_amplicon_evaluation():
    analyzer = AmpliconContextAnalyzer(target_region="16s_v3v4")

    # Mock statistics where duplication rate is high (85%)
    stats = {
        "total_reads": 1000,
        "duplication_rate": 85.0,
        "mean_read_length": 300,
        "min_read_length": 290,
        "max_read_length": 310,
        "per_base_content": {
            "A": [25.0] * 300,
            "C": [25.0] * 300,
            "G": [25.0] * 300,
            "T": [25.0] * 300,
        }
    }

    # Sample sequence starting with 341F 16S primer CCTACGGGNGGCWGCAG -> CCTACGGGAGGCAGCAG
    sample_seqs = ["CCTACGGGAGGCAGCAGTAGGGAATATTGG" * 10] * 100

    result = analyzer.evaluate(stats, sample_seqs)

    assert result["context"] == "Amplicon (16s_v3v4)"
    assert result["overall_status"] in ["PASS", "WARN"]

    # Duplication check should be PASS (expected for amplicons)
    dup_check = next(c for c in result["checks"] if "Duplication" in c["name"])
    assert dup_check["status"] == "PASS"
    assert "Standard FastQC flags high duplication" in dup_check["context_reasoning"]


def test_amplicon_long_reads_evaluation():
    analyzer = AmpliconContextAnalyzer(target_region="auto", is_long_read=True)

    stats = {
        "total_reads": 500,
        "duplication_rate": 70.0,
        "mean_read_length": 1450,
        "min_read_length": 1200,
        "max_read_length": 1650,
        "per_base_content": {
            "A": [25.0] * 1450,
            "C": [25.0] * 1450,
            "G": [25.0] * 1450,
            "T": [25.0] * 1450,
        }
    }

    # Full-length 16S Nanopore/PacBio read starting with 27F (AGAGTTTGATCMTGGCTCAG)
    sample_seqs = ["AGAGTTTGATCATGGCTCAGAGACGAACGCTGGCGGCAGGC" + "A" * 1400] * 50

    result = analyzer.evaluate(stats, sample_seqs)

    assert "Amplicon" in result["context"]
    assert result["metadata"]["sequencing_technology"].startswith("Long-Read")
    assert "16S rRNA" in result["metadata"]["target_region"]
    assert len(result["metadata"]["detected_primers"]) > 0
    assert "27F" in result["metadata"]["detected_primers"][0]["primer_name"]
