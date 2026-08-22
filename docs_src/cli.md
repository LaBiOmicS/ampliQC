# CLI Reference 💻

Complete reference for **ampliQC** command line interface.

---

## `ampliqc` Base Command

```text
Usage: ampliqc [OPTIONS] COMMAND [ARGS]...

  ampliQC: Context-aware Quality Control engine for Amplicon sequencing data (Short & Long Reads).

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  run        Analyzes FASTQ files (Single/Paired-End, Short/Long Reads).
  run-batch  Multi-Sample Parallel Batch processing across directory.
```

---

## `ampliqc run`

Analyzes single or paired-end FASTQ files.

### Syntax

```bash
ampliqc run [OPTIONS]
```

### Options

| Option | Short | Type | Default | Description |
| :--- | :--- | :--- | :--- | :--- |
| `--read1` | `-1` | Path | *Required* | Input FASTQ file (R1, single-end, or long-read). |
| `--read2` | `-2` | Path | Optional | Optional R2 FASTQ file for paired-end. |
| `--outdir` | `-o` | Path | `ampliqc_results` | Output directory for reports. |
| `--context` | `-c` | Choice | `auto` | Profiling mode: `auto`, `amplicon`, `amplicon_16s`, `amplicon_18s`, `amplicon_its`, `amplicon_long`, `ont`, `pacbio_hifi`. |
| `--primer-f` | | String | None | Custom forward primer sequence (IUPAC supported). |
| `--primer-r` | | String | None | Custom reverse primer sequence (IUPAC supported). |
| `--primer-db` | | Path | None | Custom primer database YAML/JSON file. |
| `--sample-name` | `-s` | String | Auto | Sample display name. |

---

## `ampliqc run-batch`

Runs multi-sample parallel batch processing across a directory.

### Syntax

```bash
ampliqc run-batch [OPTIONS]
```

### Options

| Option | Short | Type | Default | Description |
| :--- | :--- | :--- | :--- | :--- |
| `--input-dir` | `-i` | Path | *Required* | Input directory containing FASTQ files. |
| `--outdir` | `-o` | Path | `ampliqc_batch_results` | Output directory for batch reports. |
| `--threads` | `-t` | Int | `4` | Number of parallel worker threads. |
| `--context` | `-c` | Choice | `auto` | Amplicon context mode profile. |
| `--primer-f` | | String | None | Custom forward primer sequence. |
| `--primer-r` | | String | None | Custom reverse primer sequence. |
| `--primer-db` | | Path | None | Custom primer database YAML/JSON file. |
