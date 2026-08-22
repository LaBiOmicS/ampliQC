---
title: 'ampliQC: A Context-Aware Quality Control Engine for Short and Long-Read Amplicon Sequencing'
tags:
  - amplicon
  - quality control
  - 16S rRNA
  - 18S rRNA
  - ITS
  - metabarcoding
  - DADA2
  - Oxford Nanopore
  - PacBio HiFi
authors:
  - name: Fabiano Menegidio
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
 - name: LaBiOmicS - Laboratório de Bioinformática e Ômicas, Federal University of ABC (UFABC), Brazil
   index: 1
date: 22 August 2026
bibliography: paper.bib
---

# Summary

High-throughput amplicon sequencing (16S, 18S, ITS, and marker gene metabarcoding) is a cornerstone of microbial ecology, environmental eDNA, and clinical microbiome research. Standard Quality Control (QC) tools such as FastQC and Falco evaluate sequence reads against whole-genome sequencing assumptions, frequently flagging expected amplicon properties—such as sequence duplication and initial base composition bias—as quality failures. **ampliQC** is an open-source, context-aware QC engine written in Python that evaluates amplicon library quality within its biological and technical context for both short-read (Illumina) and full-length long-read (PacBio HiFi, Oxford Nanopore) platforms.

# Statement of Need

Generic QC engines do not account for PCR primer annealing signatures, high locus duplication rates, or amplicon target length bounds. **ampliQC** addresses this gap by integrating:

1. **Amplicon-Aware Biological Rules**: Distinguishes expected PCR amplification duplication from library bottlenecks.
2. **Universal Primer & Adapter Scanner**: Scans 5' and 3' regions using IUPAC degenerate codes and orientation detection.
3. **Automated Cutadapt & DADA2 Parameter Recommendation**: Calculates optimal Phred truncation points and generates copyable CLI commands for downstream ASV clustering.
4. **Multi-Sample Parallel Batch Processing & MultiQC Export**: Supports `ampliqc run-batch` and generates native MultiQC JSON custom content files (`_ampliqc_mqc.json`).

# References
