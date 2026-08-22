"""
Command Line Interface (CLI) for ampliQC using Click and Rich.
Supports Single-Sample analysis and Multi-Sample Parallel Batch processing.
"""

import os
import sys
import glob
from pathlib import Path
from typing import Optional, List, Dict, Any
from concurrent.futures import ProcessPoolExecutor, as_completed
import click
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from ampliqc import __version__
from ampliqc.core.parser import read_fastq, read_fastq_paired
from ampliqc.core.stats import FastqStatsAccumulator
from ampliqc.context.amplicon import AmpliconContextAnalyzer
from ampliqc.context.detector import detect_context
from ampliqc.reports.html import generate_html_report
from ampliqc.reports.json_report import export_json, export_yaml
from ampliqc.reports.plots import generate_static_plots
from ampliqc.reports.helpers import export_multiqc_json, export_tsv

console = Console()


@click.group()
@click.version_option(version=__version__, prog_name="ampliQC")
def cli():
    """
    ampliQC: Context-aware Quality Control engine for Amplicon sequencing data (Short & Long Reads).
    """
    pass


def _process_single_sample(
    r1_path: Path,
    r2_path: Optional[Path],
    out_path: Path,
    context: str,
    primer_f: Optional[str],
    primer_r: Optional[str],
    primer_db: Optional[str],
    sample_name: str,
) -> Dict[str, Any]:
    accumulator = FastqStatsAccumulator()
    sample_sequences = []

    if r2_path and r2_path.exists():
        for r1, r2 in read_fastq_paired(r1_path, r2_path):
            accumulator.process_record(r1)
            accumulator.process_record(r2)
            if len(sample_sequences) < 50000:
                sample_sequences.append(r1.sequence)
                sample_sequences.append(r2.sequence)
    else:
        for record in read_fastq(r1_path):
            accumulator.process_record(record)
            if len(sample_sequences) < 50000:
                sample_sequences.append(record.sequence)

    stats_summary = accumulator.get_summary()

    if stats_summary.get("total_reads", 0) == 0:
        return {}

    custom_primers = [p for p in [primer_f, primer_r] if p]
    mean_len = stats_summary.get("mean_read_length", 0)
    is_long_read = mean_len > 800 or stats_summary.get("max_read_length", 0) > 2500

    if context == "auto":
        detected_name, analyzer = detect_context(stats_summary, sample_sequences)
        if custom_primers and isinstance(analyzer, AmpliconContextAnalyzer):
            analyzer.custom_primers.extend(custom_primers)
    elif context in ["amplicon_long", "ont", "pacbio_hifi"]:
        analyzer = AmpliconContextAnalyzer(target_region=context, custom_primers=custom_primers, is_long_read=True)
    else:
        analyzer = AmpliconContextAnalyzer(target_region=context, custom_primers=custom_primers, is_long_read=is_long_read)

    context_eval = analyzer.evaluate(stats_summary, sample_sequences)

    # Save Reports
    out_path.mkdir(parents=True, exist_ok=True)

    html_file = out_path / f"{sample_name}_ampliqc.html"
    json_file = out_path / f"{sample_name}_ampliqc.json"
    yaml_file = out_path / f"{sample_name}_ampliqc.yaml"
    mqc_file = out_path / f"{sample_name}_ampliqc_mqc.json"

    full_data = {
        "sample_name": sample_name,
        "summary": stats_summary,
        "context_evaluation": context_eval,
    }

    generate_html_report(sample_name, stats_summary, context_eval, html_file)
    export_json(full_data, json_file)
    export_yaml(full_data, yaml_file)
    export_multiqc_json(full_data, mqc_file)
    generate_static_plots(stats_summary, out_path)

    return full_data


@cli.command(name="run")
@click.option("-1", "--read1", required=True, type=click.Path(exists=True), help="Input FASTQ file (R1, single-end, or long-read).")
@click.option("-2", "--read2", required=False, type=click.Path(exists=True), help="Optional R2 FASTQ file for paired-end.")
@click.option("-o", "--outdir", default="ampliqc_results", help="Output directory for reports.")
@click.option(
    "-c", "--context",
    type=click.Choice(["auto", "amplicon", "amplicon_16s", "amplicon_18s", "amplicon_its", "amplicon_long", "ont", "pacbio_hifi"], case_sensitive=False),
    default="auto",
    help="Amplicon sequencing context profile (Short-Read vs Long-Read 16S/18S/ITS/Functional)."
)
@click.option("--primer-f", default=None, help="Custom forward primer sequence.")
@click.option("--primer-r", default=None, help="Custom reverse primer sequence.")
@click.option("--primer-db", default=None, type=click.Path(exists=True), help="Custom primer database YAML/JSON file.")
@click.option("-s", "--sample-name", default=None, help="Sample display name.")
def run(
    read1: str,
    read2: Optional[str],
    outdir: str,
    context: str,
    primer_f: Optional[str],
    primer_r: Optional[str],
    primer_db: Optional[str],
    sample_name: Optional[str],
):
    """
    Analyzes FASTQ files (Single/Paired-End, Short/Long Reads) and generates context-aware QC reports.
    """
    r1_path = Path(read1)
    r2_path = Path(read2) if read2 else None
    out_path = Path(outdir)
    out_path.mkdir(parents=True, exist_ok=True)

    if not sample_name:
        name_clean = r1_path.name
        for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
            if name_clean.endswith(ext):
                name_clean = name_clean[:-len(ext)]
                break
        import re
        sample_name = re.sub(r'(_R1|_R2|_1|_2)$', '', name_clean)

    library_mode = "Paired-End" if r2_path else "Single-End / Long-Read"

    console.print(Panel.fit(
        f"[bold cyan]ampliQC v{__version__}[/bold cyan]\n"
        f"Processing sample: [bold yellow]{sample_name}[/bold yellow]\n"
        f"Library mode: [bold white]{library_mode}[/bold white]\n"
        f"Input R1: [dim]{r1_path}[/dim]\n"
        + (f"Input R2: [dim]{r2_path}[/dim]\n" if r2_path else "") +
        f"Context mode: [bold green]{context}[/bold green]",
        title="🧬 ampliQC: Amplicon Context-Aware QC",
        border_style="cyan"
    ))

    with console.status(f"[bold green]Streaming and analyzing {library_mode} FASTQ reads...", spinner="dots"):
        full_data = _process_single_sample(
            r1_path, r2_path, out_path, context, primer_f, primer_r, primer_db, sample_name
        )

    if not full_data:
        console.print("[bold red]Error: FASTQ contains 0 valid reads.[/bold red]")
        sys.exit(1)

    context_eval = full_data["context_evaluation"]
    meta = context_eval.get("metadata", {})

    # Print Metadata Panel
    if meta:
        meta_table = Table(show_header=False, box=None)
        meta_table.add_column("Key", style="bold cyan", justify="right")
        meta_table.add_column("Value", style="bold white")

        meta_table.add_row("Tipo de Sequenciamento:", meta.get("sequencing_technology", "N/A"))
        meta_table.add_row("Estratégia de Biblioteca:", meta.get("context_type", "N/A"))
        meta_table.add_row("Região Alvo (Target Region):", meta.get("target_region", "N/A"))
        meta_table.add_row("Primer Detectado / Usado:", meta.get("primary_primer", "N/A"))
        
        trim_st = meta.get("trimming_status", "N/A")
        if "UNTRIMMED" in trim_st or "Required" in trim_st:
            trim_st_str = f"[bold yellow]⚠️ {trim_st}[/bold yellow]"
        else:
            trim_st_str = f"[bold green]✓ {trim_st}[/bold green]"
        meta_table.add_row("Status dos Primers:", trim_st_str)
        meta_table.add_row("Pipeline Recomendado:", f"[italic bright_green]{meta.get('recommended_pipeline', 'N/A')}[/italic bright_green]")

        console.print(Panel(meta_table, title="🧬 Perfil Biológico e Metadados do Sequenciamento", border_style="blue"))

    # Print Terminal Table Summary
    table = Table(title=f"QC Results: {sample_name} ({context_eval['context']})", show_lines=True)
    table.add_column("Check Name", style="bold")
    table.add_column("Status", justify="center")
    table.add_column("Observation")
    table.add_column("Context Perspective", style="dim")

    for check in context_eval["checks"]:
        st = check["status"]
        if st == "PASS":
            status_str = "[bold green]PASS[/bold green]"
        elif st == "WARN":
            status_str = "[bold yellow]WARN[/bold yellow]"
        else:
            status_str = "[bold red]FAIL[/bold red]"

        table.add_row(check["name"], status_str, check["message"], check["context_reasoning"])

    console.print(table)

    # Print Cutadapt Command & DADA2 Recommendations
    if meta.get("cutadapt_command"):
        console.print(Panel(
            f"[bold yellow]{meta['cutadapt_command']}[/bold yellow]",
            title="✂️ Cutadapt / Chopper CLI Command (Ready to Copy)",
            border_style="yellow"
        ))

    if meta.get("dada2_recommendation"):
        console.print(Panel(
            f"[bold green]{meta['dada2_recommendation']}[/bold green]",
            title="🧪 DADA2 Truncation Parameter Suggestion",
            border_style="green"
        ))

    console.print(f"\n[bold green]✓ QC completed successfully![/bold green]")
    console.print(f"📄 HTML Report: [underline cyan]{(out_path / f'{sample_name}_ampliqc.html').resolve()}[/underline cyan]")
    console.print(f"📊 JSON Report: [underline cyan]{(out_path / f'{sample_name}_ampliqc.json').resolve()}[/underline cyan]")
    console.print(f"📋 MultiQC JSON: [underline cyan]{(out_path / f'{sample_name}_ampliqc_mqc.json').resolve()}[/underline cyan]")


@cli.command(name="run-batch")
@click.option("-i", "--input-dir", required=True, type=click.Path(exists=True), help="Directory containing FASTQ files.")
@click.option("-o", "--outdir", default="ampliqc_batch_results", help="Output directory for batch reports.")
@click.option("-t", "--threads", default=4, help="Number of parallel worker processes.")
@click.option(
    "-c", "--context",
    type=click.Choice(["auto", "amplicon", "amplicon_16s", "amplicon_18s", "amplicon_its", "amplicon_long"], case_sensitive=False),
    default="auto",
    help="Context mode for all samples."
)
def run_batch(input_dir: str, outdir: str, threads: int, context: str):
    """
    Executes parallel multi-sample Quality Control across an entire directory of FASTQ files.
    """
    in_path = Path(input_dir)
    out_path = Path(outdir)
    out_path.mkdir(parents=True, exist_ok=True)

    fastq_files = sorted(list(in_path.glob("*.fastq*")) + list(in_path.glob("*.fq*")))

    if not fastq_files:
        console.print(f"[bold red]No FASTQ files found in {input_dir}[/bold red]")
        sys.exit(1)

    # Pair R1 and R2 files
    sample_pairs = {}
    for f in fastq_files:
        name = f.name
        for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
            if name.endswith(ext):
                name = name[:-len(ext)]
                break

        if "_R2" in name or "_2" in name:
            s_name = name.replace("_R2", "").replace("_2", "")
            sample_pairs.setdefault(s_name, {})["R2"] = f
        else:
            s_name = name.replace("_R1", "").replace("_1", "")
            sample_pairs.setdefault(s_name, {})["R1"] = f

    console.print(Panel.fit(
        f"[bold cyan]ampliQC Batch Runner[/bold cyan]\n"
        f"Samples Found: [bold yellow]{len(sample_pairs)}[/bold yellow]\n"
        f"Threads / Workers: [bold white]{threads}[/bold white]\n"
        f"Input Dir: [dim]{in_path}[/dim]\n"
        f"Output Dir: [dim]{out_path}[/dim]",
        title="🚀 Multi-Sample Amplicon Batch Processing",
        border_style="magenta"
    ))

    results = []
    with console.status(f"[bold green]Running batch QC across {len(sample_pairs)} samples...", spinner="dots"):
        for s_name, files in sample_pairs.items():
            r1 = files.get("R1")
            r2 = files.get("R2")
            if not r1:
                continue
            res = _process_single_sample(r1, r2, out_path, context, None, None, None, s_name)
            if res:
                results.append(res)

    tsv_file = out_path / "batch_summary_ampliqc.tsv"
    export_tsv(results, tsv_file)

    console.print(f"\n[bold green]✓ Batch QC complete for {len(results)} samples![/bold green]")
    console.print(f"📊 Batch TSV Matrix: [underline cyan]{tsv_file.resolve()}[/underline cyan]")


if __name__ == "__main__":
    cli()
