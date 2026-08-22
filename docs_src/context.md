# Biological Context & Rules 🧬

This document details the biological rules and evaluation criteria used by **ampliQC**.

---

## Amplicon Sequencing Context vs Whole Genome Sequencing

Standard Quality Control tools (such as FastQC and Falco) evaluate sequencing datasets against expectations derived from random Whole Genome Sequencing (WGS). 

In WGS:
- Reads start at random positions across the genome.
- Nucleotide composition ($A, C, G, T$) is roughly balanced.
- High sequence duplication indicates a PCR amplification bottleneck or optical duplicate error.

In Amplicon Libraries:
- Reads originate from specific, targeted loci (e.g., 16S V3-V4, ITS1, COI).
- Reads start at fixed PCR primer binding sites.
- Millions of reads share identical sequence fragments due to biological abundance.

---

## Supported Target Regions

`ampliQC` includes built-in knowledge for major marker gene amplicon targets:

| Marker Gene | Target Region | Typical Length | Common Primers |
| :--- | :--- | :--- | :--- |
| **16S rRNA** | V1-V3 | ~490 bp | 27F / 534R |
| **16S rRNA** | V3-V4 | ~460 bp | 341F / 805R |
| **16S rRNA** | V4 | ~290 bp | 515F / 806R |
| **16S rRNA** | Full-Length | ~1,500 bp | 27F / 1492R |
| **18S rRNA** | V4 / V9 | ~350–450 bp | TAReuk454FWD1 / TAReukREV3 |
| **ITS** | ITS1 / ITS2 | ~250–400 bp | ITS1F / ITS2 / ITS4 |
| **COI** | Folmer / Leray | ~313–658 bp | LCO1490 / HCO2198 / mlCOIintF |
| **12S rRNA** | MiFish | ~170 bp | MiFish-U-F / MiFish-U-R |

---

## Evaluation Checks & Thresholds

### 1. Primer Detection & Trimming Status

- **Method**: IUPAC degenerate regex scanning across the first and last 50 nucleotides of each read.
- **PASS**: Primers are absent or detected at $<5\%$ frequency (`TRIMMED`).
- **WARN**: Primers detected in $\ge 5\%$ of reads (`UNTRIMMED`). Generates Cutadapt / Chopper trimming command.

### 2. Sequence Duplication Rate

- **Method**: K-mer and full sequence duplication counting.
- **PASS**: Duplication rate $>50\%$ is expected and marked **PASS** for amplicon libraries, accompanied by context explanation.
- **WARN**: Duplication rate $<5\%$ on high-depth libraries may indicate non-specific amplification or host contamination.

### 3. Per-Base Quality Decay & DADA2 Truncation

- **Method**: Sliding window Phred score decay ($Q < 25$ or $Q < 20$).
- **Action**: Calculates exact truncation positions (`truncLen` R1/R2) ensuring sufficient read overlap ($>20\text{ bp}$) for paired-end merging in DADA2 or QIIME 2.

### 4. Read Length Profile

- **Method**: Histogram analysis of read length distribution.
- **PASS**: Uniform read length distribution matching target amplicon platform expectations.
