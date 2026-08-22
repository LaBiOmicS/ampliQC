import gzip
import pytest
from pathlib import Path
from click.testing import CliRunner
from ampliqc.cli import cli
from ampliqc.core.primers import load_custom_primer_db, scan_primers
from ampliqc.reports.helpers import (
    generate_cutadapt_command,
    calculate_dada2_recommendation,
    export_multiqc_json,
    export_tsv,
)


def test_cutadapt_command_generation():
    detected = [{"primer_name": "341F", "sequence": "CCTACGGGNGGCWGCAG"}]
    cmd_pe = generate_cutadapt_command("sample_1", "R1.fq.gz", "R2.fq.gz", detected, is_long_read=False)
    assert "cutadapt -g CCTACGGGNGGCWGCAG" in cmd_pe
    assert "-o sample_1_trimmed_R1.fastq.gz" in cmd_pe

    cmd_long = generate_cutadapt_command("sample_ont", "ont.fq.gz", None, detected, is_long_read=True)
    assert "--revcomp" in cmd_long


def test_dada2_recommendation():
    stats = {"mean_read_length": 250, "per_base_quality": {}}
    rec_pe = calculate_dada2_recommendation(stats, "16S rRNA (V3-V4 region)", is_paired=True)
    assert "truncLen=c(" in rec_pe["dada2_code"]
    assert rec_pe["trunc_len_f"] > 100


def test_custom_primer_db_loading(tmp_path):
    yaml_file = tmp_path / "my_primers.yaml"
    yaml_file.write_text("LabPrimerF: CCTACGGGAGGCAG\nLabPrimerR: GACTACHVGGGT\n")

    loaded = load_custom_primer_db(str(yaml_file))
    assert loaded["LabPrimerF"] == "CCTACGGGAGGCAG"
    assert loaded["LabPrimerR"] == "GACTACHVGGGT"


def test_multiqc_export_and_tsv(tmp_path):
    full_data = {
        "sample_name": "S1",
        "summary": {"total_reads": 1000, "mean_read_length": 250, "overall_gc_content": 52.0, "duplication_rate": 45.0},
        "context_evaluation": {
            "context": "Amplicon (16s)",
            "overall_status": "PASS",
            "metadata": {"target_region": "16S V4", "primary_primer": "515F", "trimming_status": "TRIMMED / CLEAN"}
        }
    }

    mqc_json = tmp_path / "S1_mqc.json"
    export_multiqc_json(full_data, mqc_json)
    assert mqc_json.exists()

    tsv_file = tmp_path / "batch.tsv"
    export_tsv([full_data], tsv_file)
    assert tsv_file.exists()
    content = tsv_file.read_text()
    assert "S1" in content
    assert "16S V4" in content


def test_cli_run_batch(tmp_path):
    # Create mock FASTQ files for 2 samples
    for sname in ["sampleA", "sampleB"]:
        r1_file = tmp_path / f"{sname}_R1.fastq"
        seq = "CCTACGGGAGGCAGCAGTAGGGAATATTGG" * 5
        qual = "I" * len(seq)
        with open(r1_file, "w") as f:
            for i in range(20):
                f.write(f"@{sname}_{i}\n{seq}\n+\n{qual}\n")

    out_dir = tmp_path / "batch_output"
    runner = CliRunner()
    res = runner.invoke(cli, [
        "run-batch",
        "-i", str(tmp_path),
        "-o", str(out_dir),
        "-t", "2"
    ])

    assert res.exit_code == 0, f"CLI Batch Output: {res.output}"
    assert (out_dir / "batch_summary_ampliqc.tsv").exists()
