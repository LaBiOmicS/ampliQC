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

## What is ampliQC?

**ampliQC** is a context-aware Quality Control (QC) engine specifically designed for **Amplicon Sequencing Libraries** (16S rRNA, 18S rRNA, ITS, and functional marker genes like COI, 12S, rbcL, matK, nifH, amoA). It supports both **Short-Read (Illumina / MGI / Element)** and **Long-Read (PacBio HiFi / Oxford Nanopore)** sequencing platforms.

### The Problem with Generic QC Tools

Standard quality control utilities such as FastQC or Falco were created under the assumption of **Whole Genome Sequencing (WGS)** or **RNA-Seq** libraries, where reads are randomly fragmented across the entire genome. 

When applied to amplicon libraries, generic QC tools trigger **false-positive warnings and failures**:

1. 🛑 **Sequence Duplication Failure**: FastQC flags high duplication rates as a library bottleneck error. In amplicon sequencing, millions of reads target the exact same genomic locus, making high sequence duplication a normal biological feature.
2. 🛑 **Per-Base Sequence Content Failure**: FastQC flags non-uniform base distribution ($A/C/G/T$) at read starts. Amplicon reads start at fixed PCR primer binding sites, exhibiting conserved sequence bias that is expected and necessary.
3. 🛑 **Ignored Untrimmed Primers**: FastQC does not detect remaining PCR primers or calculate downstream ASV/OTU impact.

### The ampliQC Solution

`ampliQC` evaluates quality metrics **in context**:
- It distinguishes expected target amplification from actual quality degradation.
- It detects degenerate PCR primers (IUPAC codes) at read ends and checks trimming status.
- It automatically generates ready-to-use **Cutadapt** / **Chopper** commands.
- It calculates optimal **DADA2 / QIIME 2 truncation lengths** (`truncLen`).

---

## Key Features

- 🎯 **Amplicon-Aware Biological Evaluation**: 16S (V1-V3, V3-V4, V4, V9, Full-Length), 18S, ITS, COI, 12S, rbcL, matK, nifH, amoA.
- ✂️ **Automated Trimming & DADA2 Recommendations**: Instant Cutadapt / Chopper CLI generation and DADA2 `truncLen` calculation.
- 🧬 **Degenerate PCR Primer & Flank Detection**: IUPAC regex scanning across 5' and 3' read ends.
- 🚀 **Multi-Sample Parallel Batch Engine**: Fast directory processing via `ampliqc run-batch` with MultiQC JSON and TSV matrix exports.
- 🛠️ **Custom Primer Panel Support**: YAML configuration support for lab-specific primer sets.
- 📊 **Rich Multi-Format Reporting**: Standalone HTML, structured JSON/YAML, MultiQC JSON, and static PNG figures.

---

## Quick Example

```bash
# Analyze a 16S amplicon dataset (Illumina V3-V4 Paired-End)
ampliqc run -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz --context auto -o qc_results/

# Run multi-sample batch QC in parallel (8 threads)
ampliqc run-batch -i /path/to/fastq_dir/ -o batch_qc/ -t 8
```
