# ampliQC 🧬

**ampliQC** is a context-aware Quality Control (QC) engine specialized for **Amplicon Sequencing Data** (16S, 18S, ITS, and functional marker genes) across both **Short-Read (Illumina / MGI / Element)** and **Long-Read (PacBio HiFi / Oxford Nanopore)** platforms.

Unlike generic QC tools like FastQC or Falco—which evaluate amplicon reads against genome sequencing assumptions and trigger false warnings for high duplication rates or initial primer nucleotide composition—**ampliQC** evaluates sequencing quality **within the biological and technical context of amplicon libraries**.

---

## Key Features

- 🎯 **Amplicon-Aware Biological Evaluation**:
  - **16S / 18S / ITS & Custom Target Profiling**: Detects target regions (V1-V3, V3-V4, V4, V9, ITS1, ITS2, or Full-Length ~1.5kb).
  - **Primer Detection & Trimming Status**: Scans 5' and 3' ends for degenerate PCR primers (IUPAC codes) and assesses trimming state.
  - **Context-Aware Duplication Rate**: Distinguishes expected biological locus amplification from library bottleneck errors.
  - **Conserved Base Bias vs Body Composition**: Separates expected primer binding site bias from target sequence quality drops.
  - **Short vs Long-Read Profiling**: Supports Illumina (150-300bp) and PacBio HiFi / ONT Full-Length amplicons.
- ✂️ **Automatic Trimming & DADA2 Parameters**:
  - Generates ready-to-copy **Cutadapt** and **Chopper** CLI commands.
  - Calculates optimal **DADA2 / QIIME 2 truncation lengths** based on Phred quality decay and amplicon target overlap.
- 🚀 **Multi-Sample Parallel Batch Processing**:
  - Process an entire directory of FASTQ files in parallel via `ampliqc run-batch`.
  - Generates consolidated Multi-Sample **TSV summary matrices** and native **MultiQC JSON** files (`_ampliqc_mqc.json`).
- 🛠️ **Custom Primer Database Support**:
  - Load custom lab primer panels via `--primer-db custom_primers.yaml` or default `~/.config/ampliqc/primers.yaml`.

---

## Installation

```bash
pip install -e .
```

## Quick Start

```bash
# Analyze a 16S amplicon dataset (Paired-End Illumina V3-V4)
ampliqc run -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz --context amplicon_16s -o qc_output/

# Analyze a Full-Length 16S dataset (PacBio HiFi / Oxford Nanopore Long-Read)
ampliqc run -1 sample_longreads.fastq.gz --context auto -o qc_output/

# Run multi-sample batch QC across a directory in parallel (4 threads)
ampliqc run-batch -i /path/to/fastq_dir/ -o batch_qc_output/ -t 4

# Specify custom forward and reverse primers or custom YAML database
ampliqc run -1 sample.fastq.gz --primer-f GTGYCAGCMGCCGCGGTAA --primer-r GGACTACNVGGGTWTCTAAT -o qc_output/
ampliqc run -1 sample.fastq.gz --primer-db lab_primers.yaml -o qc_output/
```
