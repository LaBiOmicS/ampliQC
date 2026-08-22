# ampliQC 🧬

**ampliQC** is a context-aware Quality Control (QC) engine specialized for **Amplicon Sequencing Data** (16S, 18S, ITS, and functional marker genes) across both **Short-Read (Illumina / MGI / Element)** and **Long-Read (PacBio HiFi / Oxford Nanopore)** platforms.

---

## Key Features

- 🎯 **Amplicon-Aware Biological Evaluation** (16S, 18S, ITS, COI, 12S MiFish, rbcL, matK, nifH, amoA).
- ✂️ **Automatic Trimming & DADA2 Parameters** (Cutadapt & Chopper CLI generator, DADA2 truncation length calculator).
- 🚀 **Multi-Sample Parallel Batch Processing** (`ampliqc run-batch`).
- 📊 **MultiQC & TSV Exports** (`_ampliqc_mqc.json` and TSV summary matrix).
- 🛠️ **Custom Primer Database Support** (`--primer-db custom_primers.yaml`).
