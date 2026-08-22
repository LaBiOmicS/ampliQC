import gzip
import pytest
from click.testing import CliRunner
from ampliqc.cli import cli


def test_cli_run_amplicon(tmp_path):
    # Create sample gzipped FASTQ file
    fastq_gz = tmp_path / "sample_16s_R1.fastq.gz"
    out_dir = tmp_path / "results"

    seq = "CCTACGGGAGGCAGCAGTAGGGAATATTGGACAATGGGCGCAAGCCTGATCCAGCCATGCCGCGTGTGTGAAGAAGGCCTTCGGGTTGTAAAGCACTTT"
    qual = "I" * len(seq)

    with gzip.open(fastq_gz, "wt") as f:
        for i in range(100):
            f.write(f"@read_{i}\n{seq}\n+\n{qual}\n")

    runner = CliRunner()
    res = runner.invoke(cli, [
        "run",
        "-1", str(fastq_gz),
        "-o", str(out_dir),
        "-c", "amplicon_16s"
    ])

    assert res.exit_code == 0, f"CLI output: {res.output}"
    assert (out_dir / "sample_16s_ampliqc.html").exists()
    assert (out_dir / "sample_16s_ampliqc.json").exists()
    assert (out_dir / "sample_16s_ampliqc.yaml").exists()
