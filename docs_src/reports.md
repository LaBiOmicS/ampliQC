# Reports & Deliverables 📊

This page describes all output deliverables produced by **ampliQC**.

---

## Deliverables Summary

Each execution of `ampliQC` creates a dedicated output folder:

```
qc_output/
├── sample_ampliqc.html         # Interactive standalone HTML report
├── sample_ampliqc.json         # Raw metrics & context evaluation (JSON)
├── sample_ampliqc.yaml         # Raw metrics & context evaluation (YAML)
├── sample_ampliqc_mqc.json     # Native MultiQC integration module JSON
└── images/                     # FastQC/Falco compatible PNG images
    ├── per_base_quality.png
    ├── per_base_sequence_content.png
    ├── per_sequence_quality.png
    ├── per_sequence_gc_content.png
    └── sequence_length_distribution.png
```

For batch processing (`ampliqc run-batch`), a multi-sample summary matrix is additionally exported:
- `ampliqc_batch_summary.tsv`

---

## 1. Interactive HTML Report (`*_ampliqc.html`)

The HTML report is completely standalone (CSS and JavaScript embedded via Jinja2 templates) and can be viewed directly in any web browser without internet connectivity.

Key Sections:
- **Header & Summary Badges**: Overall sample status, read counts, read length, GC content, target region.
- **Sequencing & Primer Metadata Panel**: Technology, target amplicon region, primary primer detected, trimming status.
- **Context-Aware Evaluation Table**: List of checks, PASS/WARN/FAIL status, observed values, biological reasoning.
- **Interactive Plotly Charts**: Interactive per-base quality curves, base content, sequence length distribution.
- **Cutadapt / Chopper Ready-to-Copy Panel**: Pre-formatted CLI commands with exact primer sequences.
- **DADA2 Parameter Suggestion**: Recommended `truncLen`, `maxEE`, `truncQ` settings for ASV clustering.

---

## 2. MultiQC Integration (`*_ampliqc_mqc.json`)

`ampliQC` outputs native MultiQC custom data files compatible with MultiQC reports:

```json
{
  "id": "ampliqc_module",
  "plot_type": "generalstats",
  "data": {
    "sample_1": {
      "total_reads": 50000,
      "mean_read_length": 301.0,
      "primer_trimming_status": "UNTRIMMED"
    }
  }
}
```

---

## 3. Raw Data Exports (`.json` & `.yaml`)

Structured JSON and YAML files contain all raw metrics, per-base quality distributions, k-mer counts, primer detection statistics, and evaluation check logs for downstream programmatic analysis.
