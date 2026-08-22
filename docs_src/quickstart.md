# Quick Start Guide 🚀

This guide provides step-by-step instructions for installing and running **ampliQC**.

---

## Installation

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
pip install -e .
```

### Option 4: Docker Container

```bash
docker build -t ampliqc .
docker run --rm -v $(pwd):/data ampliqc run -1 /data/sample_R1.fastq.gz -2 /data/sample_R2.fastq.gz -o /data/qc_output
```

### Option 5: Apptainer / Singularity (HPC)

```bash
apptainer build ampliqc.sif Apptainer.def
apptainer run ampliqc.sif run -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz -o qc_output/
```

---

## basic Usage

### Single-Sample Analysis (Short-Read Illumina / MGI)

Analyze a paired-end amplicon library with automatic context detection:

```bash
ampliqc run \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  --context auto \
  -o qc_output/
```

### Single-Sample Analysis (Long-Read PacBio HiFi / Nanopore)

Analyze full-length 16S / ITS long-read amplicons:

```bash
# Oxford Nanopore full-length 16S
ampliqc run -1 sample_ont.fastq.gz --context ont -o qc_output/

# PacBio HiFi full-length 16S
ampliqc run -1 sample_pacbio.fastq.gz --context pacbio_hifi -o qc_output/
```

### Custom Primers & Custom Database

Provide custom forward and reverse primer sequences directly on the command line or via a custom YAML database:

```bash
# Direct IUPAC primer input
ampliqc run -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  --primer-f GTGYCAGCMGCCGCGGTAA \
  --primer-r GGACTACNVGGGTWTCTAAT \
  -o qc_output/

# Load custom YAML primer database
ampliqc run -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  --primer-db my_primers.yaml \
  -o qc_output/
```

### Multi-Sample Parallel Batch Processing

Process an entire directory containing dozens or hundreds of FASTQ files using multi-threading:

```bash
ampliqc run-batch \
  -i /path/to/fastq_directory/ \
  -o batch_qc_results/ \
  -t 8 \
  --context auto
```
