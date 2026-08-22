# ampliQC 🧬

<p align="center">
  <strong>Context-aware Quality Control engine for Amplicon sequencing data (Short & Long Reads)</strong>
</p>

<!-- Institutional Badges -->
[![University: UMC](https://img.shields.io/badge/University-UMC-0D47A1.svg)](https://www.umc.br/)
[![Laboratory: LaBiOmicS](https://img.shields.io/badge/Laboratory-LaBiOmicS-7B1FA2.svg)](https://github.com/LaBiOmicS)
[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-AmpliconQC-green.svg)](https://github.com/LaBiOmicS/ampliQC)

<!-- Open Science Badges -->
[![Open Source](https://img.shields.io/badge/Open-Source-brightgreen.svg)](https://github.com/LaBiOmicS/ampliQC)
[![Open Science](https://img.shields.io/badge/Open-Science-blue.svg)](https://github.com/LaBiOmicS/ampliQC)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<!-- Package & Python Badges -->
[![PyPI Package](https://img.shields.io/badge/PyPI-v0.1.1-blue.svg)](https://pypi.org/project/ampliqc/)
[![Python Versions](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12%20%7C%203.13-blue.svg)](https://pypi.org/project/ampliqc/)
[![Conda Package](https://img.shields.io/badge/bioconda-v0.1.1-green.svg)](https://anaconda.org/bioconda/ampliqc)

---

**ampliQC** is a next-generation, context-aware Quality Control (QC) engine specialized for **Amplicon Sequencing Data** (16S rRNA, 18S rRNA, ITS, and functional marker genes) across both **Short-Read (Illumina / MGI / Element)** and **Long-Read (PacBio HiFi / Oxford Nanopore)** platforms.

Unlike generic QC tools (e.g. FastQC, Falco)—which evaluate amplicon libraries against random genomic sequencing assumptions and trigger false warnings for high duplication rates or initial primer nucleotide composition—**ampliQC** evaluates sequencing quality **within the biological and technical context of amplicon libraries**.

---

## 🌟 Key Features

- 🎯 **Amplicon-Aware Biological QC Engine**:
  - **16S / 18S / ITS & Custom Target Profiling**: Detects hypervariable target regions (V1-V3, V3-V4, V4, V9, ITS1, ITS2, or Full-Length ~1.5kb).
  - **Context-Aware Duplication Rate**: Distinguishes expected biological locus amplification from library bottleneck errors.
  - **Conserved Base Bias vs Body Composition**: Separates expected primer binding site bias from target sequence quality drops.
  - **Short vs Long-Read Profiling**: Supports Illumina (150-300bp) and PacBio HiFi / ONT Full-Length amplicons.
- 🧬 **Degenerate PCR Primer & Flank Detection**:
  - Scans 5' and 3' ends for degenerate PCR primers using full IUPAC ambiguous nucleotide regex matching.
  - Evaluates primer trimming status (`TRIMMED` vs `UNTRIMMED`).
- ✂️ **Automatic Trimming & DADA2 Parameters**:
  - Generates ready-to-copy **Cutadapt** and **Chopper** CLI commands with exact primer sequences.
  - Calculates optimal **DADA2 / QIIME 2 truncation lengths** (`truncLen` R1/R2) based on Phred quality decay and amplicon target overlap.
- 🚀 **Multi-Sample Parallel Batch Processing**:
  - Process an entire directory of FASTQ files in parallel with `ampliqc run-batch`.
  - Generates consolidated Multi-Sample **TSV summary matrices** and native **MultiQC JSON** files (`_ampliqc_mqc.json`).
- 🛠️ **Custom Primer Database Support**:
  - Load custom lab primer panels via `--primer-db custom_primers.yaml` or default `~/.config/ampliqc/primers.yaml`.
- 🔬 **FAIR Compliance & Containerized Reproducibility**:
  - Includes Dockerfile and Apptainer/Singularity manifests for 100% deterministic execution in HPC and Cloud environments.

---

## 📦 Installation

### Option 1: Via PyPI (Recommended)

```bash
pip install ampliqc
```

### Option 2: Via Conda / Bioconda

```bash
conda install -c bioconda ampliqc
# or using mamba
mamba install -c bioconda ampliqc
```

### Option 3: From Source (GitHub)

```bash
git clone https://github.com/LaBiOmicS/ampliQC.git
cd ampliQC

# Install in editable mode
pip install -e .
```

### Option 4: Docker Container

```bash
docker build -t ampliqc .
docker run --rm -v $(pwd):/data ampliqc run -1 /data/sample_R1.fastq.gz -2 /data/sample_R2.fastq.gz -o /data/qc_output
```

### Option 5: Apptainer / Singularity (HPC Environments)

```bash
apptainer build ampliqc.sif Apptainer.def
apptainer run ampliqc.sif run -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -o qc_output/
```

---

## 🚀 Quick Start & CLI Usage

### 1. Single-Sample Analysis (Paired-End Illumina V3-V4)

```bash
# Automatic context detection (Short-Read 16S V3-V4)
ampliqc run -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -o qc_output/

# Explicit context setting
ampliqc run -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz --context amplicon_16s -o qc_output/
```

### 2. Full-Length Amplicon Analysis (PacBio HiFi / Oxford Nanopore)

```bash
# Long-Read Full-Length 16S / ITS amplicon analysis
ampliqc run -1 sample_ont_16s.fastq.gz --context ont -o qc_output/

# PacBio HiFi Full-Length 16S
ampliqc run -1 sample_pacbio.fastq.gz --context pacbio_hifi -o qc_output/
```

### 3. Custom Primers & Custom Database

```bash
# Specify custom forward and reverse primer sequences (IUPAC codes supported)
ampliqc run -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  --primer-f GTGYCAGCMGCCGCGGTAA \
  --primer-r GGACTACNVGGGTWTCTAAT \
  -o qc_output/

# Use a custom YAML primer database file
ampliqc run -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  --primer-db my_lab_primers.yaml \
  -o qc_output/
```

### 4. Multi-Sample Parallel Batch Processing

```bash
# Run multi-sample batch QC across a directory using 8 parallel CPU workers
ampliqc run-batch -i /path/to/fastq_dir/ -o batch_qc_results/ -t 8 --context auto
```

---

## 📂 Output Deliverables

Upon completion, `ampliQC` generates a structured output folder containing interactive reports, static figures, and structured data files:

```
qc_output/
├── sample_ampliqc.html         # Interactive standalone HTML report
├── sample_ampliqc.json         # Complete raw QC metrics & context evaluation (JSON)
├── sample_ampliqc.yaml         # Complete raw QC metrics & context evaluation (YAML)
├── sample_ampliqc_mqc.json     # Native MultiQC integration module JSON
└── images/                     # FastQC/Falco-compatible static PNG figures
    ├── per_base_quality.png
    ├── per_base_sequence_content.png
    ├── per_sequence_quality.png
    ├── per_sequence_gc_content.png
    └── sequence_length_distribution.png
```

For batch runs (`run-batch`), an additional consolidated matrix `ampliqc_batch_summary.tsv` is created in the output directory.

---

## 📂 Project Structure

```
ampliQC/
├── pyproject.toml               # Package configuration & dependencies (PEP 621)
├── MANIFEST.in                  # Package distribution manifest
├── README.md                    # Package documentation
├── LICENSE                      # MIT License
├── Dockerfile                   # Docker container manifest
├── Apptainer.def                # HPC Singularity / Apptainer manifest
├── .github/                     # GitHub Workflows
│   └── workflows/
│       ├── ci.yml               # Automated Pytest CI workflow
│       └── pypi.yml             # PyPI Release & Publish workflow
├── src/                         # Package source code
│   └── ampliqc/
│       ├── __init__.py
│       ├── cli.py               # Click & Rich CLI entrypoints
│       ├── context/             # Context-aware evaluation engines
│       │   ├── base.py
│       │   ├── amplicon.py
│       │   ├── longreads.py
│       │   └── detector.py
│       ├── core/                # FASTQ parsing, k-mer & primer stats
│       │   ├── parser.py
│       │   ├── stats.py
│       │   ├── primers.py
│       │   ├── adapters.py
│       │   └── kmer.py
│       └── reports/             # HTML, JSON, MultiQC & Matplotlib plots
│           ├── html.py
│           ├── json_report.py
│           ├── plots.py
│           ├── helpers.py
│           └── templates/
│               └── report.html.j2
└── tests/                       # Pytest test suite
    ├── test_amplicon.py
    ├── test_cli.py
    ├── test_longreads.py
    ├── test_primers.py
    └── test_stats.py
```

---

## 🧪 Testing

Run the full test suite locally with `pytest`:

```bash
pytest
```

---

## 📄 License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

---

## 🤝 Contributing & Support

Contributions are welcome! Please feel free to open an issue or submit a Pull Request on GitHub:
- **Repository**: [https://github.com/LaBiOmicS/ampliQC](https://github.com/LaBiOmicS/ampliQC)
- **Issues**: [https://github.com/LaBiOmicS/ampliQC/issues](https://github.com/LaBiOmicS/ampliQC/issues)
- **PyPI**: [https://pypi.org/project/ampliqc/](https://pypi.org/project/ampliqc/)
