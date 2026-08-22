# Python API Reference 🐍

Welcome to the Python API documentation for **ampliQC**.

---

## Core API Modules

- [**AmpliconContextAnalyzer**](context.md): Context-aware biological evaluation engine.
- [**FastqStatsAccumulator**](stats.md): Streaming FASTQ quality and base statistics engine.
- [**PrimerScanner**](primers.md): IUPAC regex primer detection and reverse complement utilities.
- [**Exporters & Reports**](reports.md): HTML, JSON, MultiQC, and Matplotlib figure generators.

---

## Simple API Example

```python
from ampliqc.core.parser import read_fastq
from ampliqc.core.stats import FastqStatsAccumulator
from ampliqc.context.amplicon import AmpliconContextAnalyzer

# 1. Accumulate FASTQ statistics
accumulator = FastqStatsAccumulator()
sequences = []

for record in read_fastq("sample.fastq.gz"):
    accumulator.process_record(record)
    if len(sequences) < 50000:
        sequences.append(record.sequence)

summary = accumulator.get_summary()

# 2. Evaluate context
analyzer = AmpliconContextAnalyzer(target_region="amplicon_16s")
evaluation = analyzer.evaluate(summary, sample_sequences=sequences)

print(f"Context: {evaluation['context']}")
for check in evaluation['checks']:
    print(f"- {check['name']}: {check['status']} ({check['message']})")
```
