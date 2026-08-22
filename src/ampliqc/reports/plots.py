"""
Static plot generator module for FastQC/Falco compatible PNG/SVG image export.
"""

from pathlib import Path
from typing import Dict, Any, Union
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def generate_static_plots(summary: Dict[str, Any], output_dir: Union[str, Path]) -> Path:
    """
    Generates static PNG plots compatible with FastQC / Falco image exports.
    Saves images into output_dir/images/.
    """
    out_images_dir = Path(output_dir) / "images"
    out_images_dir.mkdir(parents=True, exist_ok=True)

    _plot_per_base_quality(summary, out_images_dir / "per_base_quality.png")
    _plot_per_base_sequence_content(summary, out_images_dir / "per_base_sequence_content.png")
    _plot_per_sequence_quality(summary, out_images_dir / "per_sequence_quality.png")
    _plot_per_sequence_gc_content(summary, out_images_dir / "per_sequence_gc_content.png")
    _plot_sequence_length_distribution(summary, out_images_dir / "sequence_length_distribution.png")

    return out_images_dir


def _plot_per_base_quality(summary: Dict[str, Any], out_path: Path):
    means = summary.get("per_base_qual_mean", [])
    medians = summary.get("per_base_qual_median", [])
    if not means:
        return

    positions = np.arange(1, len(means) + 1)

    fig, ax = plt.subplots(figsize=(10, 5), dpi=150)
    
    # FastQC Quality Background Zones
    ax.axhspan(28, 42, color='#e6ffe6', alpha=0.9, zorder=0)  # Green (>28)
    ax.axhspan(20, 28, color='#ffffcc', alpha=0.9, zorder=0)  # Yellow (20-28)
    ax.axhspan(0, 20, color='#ffe6e6', alpha=0.9, zorder=0)   # Red (<20)

    ax.plot(positions, means, label="Mean Phred Score", color="#1d4ed8", linewidth=2.0, zorder=3)
    if medians:
        ax.plot(positions, medians, label="Median Phred Score", color="#047857", linestyle="--", linewidth=1.5, zorder=3)

    ax.set_title("Quality scores across all bases (Sanger / Illumina 1.9 encoding)", fontsize=12, fontweight='bold', pad=12)
    ax.set_xlabel("Position in read (bp)", fontsize=10)
    ax.set_ylabel("Quality Score (Phred)", fontsize=10)
    ax.set_ylim(0, 42)
    ax.set_xlim(1, max(positions))
    ax.grid(True, linestyle=':', alpha=0.6, zorder=1)
    ax.legend(loc="lower left", frameon=True, facecolor='white', framealpha=0.9)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _plot_per_base_sequence_content(summary: Dict[str, Any], out_path: Path):
    content = summary.get("per_base_content", {})
    if not content or "A" not in content:
        return

    a_pct = content["A"]
    c_pct = content["C"]
    g_pct = content["G"]
    t_pct = content["T"]
    positions = np.arange(1, len(a_pct) + 1)

    fig, ax = plt.subplots(figsize=(10, 5), dpi=150)
    ax.plot(positions, a_pct, label="% A", color="#ef4444", linewidth=1.8)
    ax.plot(positions, c_pct, label="% C", color="#3b82f6", linewidth=1.8)
    ax.plot(positions, g_pct, label="% G", color="#10b981", linewidth=1.8)
    ax.plot(positions, t_pct, label="% T", color="#f59e0b", linewidth=1.8)

    ax.set_title("Sequence content across all bases", fontsize=12, fontweight='bold', pad=12)
    ax.set_xlabel("Position in read (bp)", fontsize=10)
    ax.set_ylabel("Nucleotide Content (%)", fontsize=10)
    ax.set_ylim(0, 100)
    ax.set_xlim(1, max(positions))
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(loc="upper right", frameon=True, facecolor='white', framealpha=0.9)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _plot_per_sequence_quality(summary: Dict[str, Any], out_path: Path):
    hist = summary.get("read_qual_histogram", [])
    if not hist:
        return

    qual_scores = np.arange(len(hist))

    fig, ax = plt.subplots(figsize=(9, 4.5), dpi=150)
    ax.bar(qual_scores, hist, color="#6366f1", width=0.8, edgecolor="#4338ca")

    ax.set_title("Quality score distribution over all sequences", fontsize=12, fontweight='bold', pad=12)
    ax.set_xlabel("Mean Sequence Quality (Phred Score)", fontsize=10)
    ax.set_ylabel("Number of Reads", fontsize=10)
    ax.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _plot_per_sequence_gc_content(summary: Dict[str, Any], out_path: Path):
    gc_hist = summary.get("gc_content_histogram", [])
    if not gc_hist:
        return

    gc_pcts = np.arange(len(gc_hist))

    fig, ax = plt.subplots(figsize=(9, 4.5), dpi=150)
    ax.bar(gc_pcts, gc_hist, color="#0d9488", width=1.0, edgecolor="#0f766e")

    ax.set_title("GC distribution over all sequences", fontsize=12, fontweight='bold', pad=12)
    ax.set_xlabel("Mean GC Count (%)", fontsize=10)
    ax.set_ylabel("Number of Reads", fontsize=10)
    ax.set_xlim(0, 100)
    ax.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def _plot_sequence_length_distribution(summary: Dict[str, Any], out_path: Path):
    min_len = summary.get("min_read_length", 0)
    max_len = summary.get("max_read_length", 0)
    mean_len = summary.get("mean_read_length", 0)

    fig, ax = plt.subplots(figsize=(8, 4), dpi=150)
    if min_len == max_len:
        ax.bar([min_len], [summary.get("total_reads", 1)], color="#8b5cf6", width=2.0)
    else:
        ax.bar([min_len, mean_len, max_len], [summary.get("total_reads", 1)] * 3, color="#8b5cf6", width=5.0)

    ax.set_title("Sequence Length Distribution", fontsize=12, fontweight='bold', pad=12)
    ax.set_xlabel("Sequence Length (bp)", fontsize=10)
    ax.set_ylabel("Count", fontsize=10)
    ax.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
