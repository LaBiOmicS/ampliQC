"""
Core Primer & Adapter Database, Custom DB Loader, and Advanced Detection Engine.
Supports Bacterial 16S, Archaeal 16S, Fungal ITS, Eukaryotic 18S,
Functional Marker Genes (COI, nifH, amoA, 12S, rbcL, matK),
Technology-Specific Flanking Adapters (Illumina, ONT, PacBio HiFi/Kinnex),
Custom YAML/JSON primer database loading, and Mismatch Tolerance / Dimer detection.
"""

import os
import re
from pathlib import Path
from typing import Dict, Any, List, Optional
import yaml
import json

# IUPAC Degenerate Nucleotide Code mapping
IUPAC_DICT = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"
}

# Comprehensive Primer & Flank Database across all NGS & TGS Technologies
KNOWN_PRIMERS: Dict[str, str] = {
    # --- Bacterial 16S rRNA Primers ---
    "27F (16S V1-V2/V1-V3)": "AGAGTTTGATCMTGGCTCAG",
    "27F-PacBio (16S Full-Length)": "AGRGTTYGATYMTGGCTCAG",
    "338R (16S V1-V2)": "TGCTGCCTCCCGTAGGAGT",
    "534R (16S V1-V3)": "ATTACCGCGGCTGCTGG",
    "341F (16S V3-V4)": "CCTACGGGNGGCWGCAG",
    "515F (16S V4)": "GTGYCAGCMGCCGCGGTAA",
    "515F-Y (16S V4 Parada)": "GTGYCAGCMGCCGCGGTAA",
    "805R (16S V3-V4)": "GACTACHVGGGTATCTAATCC",
    "806R (16S V4 Apprill)": "GGACTACNVGGGTWTCTAAT",
    "926F (16S V6-V8)": "AAACTYAAAKGAATTGACGG",
    "1392R (16S V6-V8)": "ACGGGCGGTGTGTRC",
    "1492R (16S Full-Length)": "TACGGYTACCTTGTTACGACTT",

    # --- Archaeal 16S rRNA Primers ---
    "Arch344F (Archaeal 16S)": "ACGGGGYGCAGCAGGCGCGA",
    "Arch519F (Archaeal 16S V3-V4)": "CAGCMGCCGCGGTAA",
    "Arch806R (Archaeal 16S V4)": "GGACTACVSGGGTATCTAAT",
    "Arch915R (Archaeal 16S)": "GTGCTCCCCCGCCAATTCCT",
    "Arch1044R (Archaeal 16S)": "GGCCATGCACCWCCTCTC",

    # --- Fungal ITS Primers (ITS1, ITS2, Full-Length ITS) ---
    "ITS1F (ITS1 Fungi)": "CTTGGTCATTTAGAGGAAGTAA",
    "ITS1 (ITS1 Fungi)": "TCCGTAGGTGAACCTGCGG",
    "rITS1 (ITS1 Fungi)": "ACCTGCGGAAGGATCATT",
    "fITS7 (ITS2 Fungi)": "GTGARTCATCGAATCTTTG",
    "ITS2 (ITS1/ITS2)": "GCTGCGTTCTTCATCGATGC",
    "ITS3 (ITS2 Fungi)": "GCATCGATGAAGAACGCAGC",
    "ITS4 (ITS2/Full-ITS)": "TCCTCCGCTTATTGATATGC",
    "ITS4-Fun (ITS2 Fungi)": "AGCCTCCGCTTATTGATATG",

    # --- Eukaryotic 18S rRNA Primers ---
    "TAReuk454FWD1 (18S V4)": "CCAGCASCYGCGGTAATTCC",
    "TAReukREV3 (18S V4)": "ACTTTCGTTCTTGATYRA",
    "1389F (18S V9)": "TTGTACACACCGCCC",
    "1510R (18S V9)": "CCTTCYGCAGGTTCACCTAC",
    "Euk1391f (18S V9)": "GTACACACCGCCCGTC",
    "EukA (18S Full-Length)": "AACCTGGTTGATCCTGCCAGT",
    "EukB (18S Full-Length)": "TGATCCTTCTGCAGGTTCACCTAC",

    # --- Environmental Metabarcoding & Functional Genes ---
    "COI-LCO1490 (Barcoding)": "GGTCAACAAATCATAAAGATATTGG",
    "COI-HCO2198 (Barcoding)": "TAAACTTCAGGGTGACCAAAAAATCA",
    "COI-mlCOIintF (Eukaryotes)": "GGWACWGGWTGAACWGTWTAYCCYCC",
    "MiFish-U-F (12S Vertebrate/Fish)": "GTCGGTAAAACTCGTGCCAGC",
    "MiFish-U-R (12S Vertebrate/Fish)": "CATAGTGGGGTATCTAATCCCAGTTTG",
    "Mammal-12S-F (Mammalian Barcode)": "CGAGAAGACCCTATGGAGCT",
    "Mammal-12S-R (Mammalian Barcode)": "CCGAGGTCGCAAACACCTT",
    "rbcL-a_f (Plant Barcoding)": "ATGTCACCACAAACAGAGACTAAAGC",
    "matK-3F_KIM (Plant Barcoding)": "CGTACAGTACTTTTGTGTTTACGAG",
    "nifH-PolF (Nitrogenase)": "TGCGAYCCSAARGCBGACTC",
    "nifH-PolR (Nitrogenase)": "ATSGCCATCATYTCRCCGGA",
    "amoA-1F (Ammonia Oxidation)": "GGGGTTTCTACTGGTGGT",
    "amoA-2R (Ammonia Oxidation)": "CCCCTCKGSAAAGCCTTCTTC",

    # --- Technology Adapters, Overhangs & Flank Sequences ---
    "Illumina-Nextera-Fwd-Overhang": "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG",
    "Illumina-Nextera-Rev-Overhang": "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
    "Illumina-TruSeq-Universal": "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "ONT-16S-Flank-Fwd": "TTTCTGTTGGTGCTGATATTGC",
    "ONT-16S-Flank-Rev": "ACTTGCCTGTCGCTCTATACTC",
    "ONT-Adapter-Y": "GGCGTCTGCTTGGGTGTTTAACCT",
    "PacBio-16S-Fwd": "AGRGTTYGATYMTGGCTCAG",
    "PacBio-16S-Rev": "RGYTACCTTGTTACGACTT",
    "PacBio-Kinnex-Adapter": "CTACACGACGCTCTTCCGATCT",
}

# Primer to Region mapping lookup
PRIMER_REGION_MAP: Dict[str, str] = {
    "27F": "16S rRNA (V1-V2 / V1-V3 region)",
    "338R": "16S rRNA (V1-V2 region)",
    "534R": "16S rRNA (V1-V3 region)",
    "341F": "16S rRNA (V3-V4 region)",
    "805R": "16S rRNA (V3-V4 region)",
    "515F": "16S rRNA (V4 region)",
    "806R": "16S rRNA (V4 region)",
    "926F": "16S rRNA (V6-V8 region)",
    "1392R": "16S rRNA (V6-V8 region)",
    "1492R": "16S rRNA (Full-Length ~1.5kb)",
    "Arch344F": "Archaeal 16S rRNA",
    "Arch519F": "Archaeal 16S rRNA (V3-V4 region)",
    "Arch806R": "Archaeal 16S rRNA (V4 region)",
    "Arch915R": "Archaeal 16S rRNA",
    "Arch1044R": "Archaeal 16S rRNA",
    "ITS1F": "ITS1 (Fungi)",
    "ITS1": "ITS1 (Fungi)",
    "rITS1": "ITS1 (Fungi)",
    "fITS7": "ITS2 (Fungi)",
    "ITS2": "ITS1/ITS2 (Fungi)",
    "ITS3": "ITS2 (Fungi)",
    "ITS4": "ITS2 / Full-Length ITS (Fungi)",
    "TAReuk454FWD1": "18S rRNA (V4 region)",
    "TAReukREV3": "18S rRNA (V4 region)",
    "1389F": "18S rRNA (V9 region - Eukaryotes)",
    "1510R": "18S rRNA (V9 region - Eukaryotes)",
    "Euk1391f": "18S rRNA (V9 region - Eukaryotes)",
    "EukA": "18S rRNA (Full-Length)",
    "EukB": "18S rRNA (Full-Length)",
    "COI": "COI (Cytochrome c Oxidase I / Barcoding)",
    "MiFish": "12S rRNA (Fish / Vertebrate eDNA)",
    "Mammal-12S": "12S rRNA (Mammalian Barcode)",
    "rbcL": "rbcL (Plant Plastid Gene)",
    "matK": "matK (Plant Chloroplast Gene)",
    "nifH": "nifH (Nitrogen Fixation Gene)",
    "amoA": "amoA (Ammonia Oxidation Gene)",
    "Illumina-Nextera": "Illumina Amplicon Adapter Overhang",
    "Illumina-TruSeq": "Illumina TruSeq Adapter",
    "ONT-16S-Flank": "Oxford Nanopore Amplicon Flanking Adapter",
    "ONT-Adapter": "Oxford Nanopore Sequencing Adapter",
    "PacBio-16S": "PacBio HiFi / Kinnex 16S Adapter",
    "PacBio-Kinnex": "PacBio Kinnex MAS-Seq Adapter",
}


def iupac_to_regex(primer_seq: str) -> str:
    """Converts a primer sequence with IUPAC codes to regex pattern."""
    pattern = ""
    for char in primer_seq.upper():
        pattern += IUPAC_DICT.get(char, char)
    return pattern


def reverse_complement(seq: str) -> str:
    """Computes the reverse complement of a nucleotide sequence."""
    trans = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(trans)[::-1]


def load_custom_primer_db(db_path: str) -> Dict[str, str]:
    """
    Loads custom primer database from a YAML or JSON file.
    Example YAML format:
      MyPrimerF: "CCTACGGGNGGCWGCAG"
      MyPrimerR: "GACTACHVGGGTATCTAATCC"
    """
    p = Path(db_path)
    if not p.exists():
        return {}

    with open(p, "r", encoding="utf-8") as f:
        if p.suffix in [".yaml", ".yml"]:
            data = yaml.safe_load(f)
        elif p.suffix == ".json":
            data = json.load(f)
        else:
            data = {}

    if isinstance(data, dict):
        return {str(k): str(v).upper() for k, v in data.items()}
    return {}


def scan_primers(
    sequences: List[str],
    custom_primers: Optional[List[str]] = None,
    custom_db_path: Optional[str] = None,
    search_window: int = 80,
    max_offset: int = 25,
    min_pct: float = 1.0,
) -> List[Dict[str, Any]]:
    """
    Scans sequence sample for all known and custom primers/adapters.
    Inspects 5' head and 3' tail within search_window (bp) requiring positional terminal anchoring
    (start position <= max_offset at 5' end or end position <= max_offset from 3' end).
    This strictly isolates untrimmed PCR primers attached at read extremities and eliminates false positive internal hits.
    """
    if not sequences:
        return []

    n_samples = len(sequences)
    primers_to_test = dict(KNOWN_PRIMERS)

    # Load custom primer file if specified
    if custom_db_path:
        custom_from_file = load_custom_primer_db(custom_db_path)
        primers_to_test.update(custom_from_file)

    # User home config default primer db (~/.config/ampliqc/primers.yaml)
    default_cfg = Path.home() / ".config" / "ampliqc" / "primers.yaml"
    if default_cfg.exists():
        primers_to_test.update(load_custom_primer_db(str(default_cfg)))

    if custom_primers:
        for idx, cp in enumerate(custom_primers, 1):
            primers_to_test[f"Custom_Primer_{idx}"] = cp

    detected_primers = []

    for name, seq in primers_to_test.items():
        pattern = re.compile(iupac_to_regex(seq))
        rc_pattern = re.compile(iupac_to_regex(reverse_complement(seq)))

        fwd_matches = 0
        rev_matches = 0

        for s in sequences:
            head = s[:search_window]
            tail = s[-search_window:] if len(s) >= search_window else s

            m_head = pattern.search(head)
            m_rc_head = rc_pattern.search(head)
            m_tail = pattern.search(tail)
            m_rc_tail = rc_pattern.search(tail)

            # Check 5' head terminal match (positional anchor check <= max_offset)
            if m_head and m_head.start() <= max_offset:
                fwd_matches += 1
            elif m_rc_head and m_rc_head.start() <= max_offset:
                rev_matches += 1
            # Check 3' tail terminal match (positional anchor check <= max_offset)
            elif m_tail and (len(tail) - m_tail.end()) <= max_offset:
                fwd_matches += 1
            elif m_rc_tail and (len(tail) - m_rc_tail.end()) <= max_offset:
                rev_matches += 1

        total_matches = fwd_matches + rev_matches
        match_pct = round((total_matches / n_samples) * 100, 1)

        if match_pct >= min_pct:
            fwd_pct = round((fwd_matches / max(1, total_matches)) * 100, 1)
            rev_pct = round((rev_matches / max(1, total_matches)) * 100, 1)
            detected_primers.append({
                "primer_name": name,
                "sequence": seq,
                "match_pct": match_pct,
                "fwd_matches": fwd_matches,
                "rev_matches": rev_matches,
                "orientation_str": f"{fwd_pct}% Forward / {rev_pct}% Reverse",
            })

    # 1. Deduplicate exact identical sequences (e.g. 515F vs 515F-Y)
    unique_by_seq: Dict[str, Dict[str, Any]] = {}
    for p in detected_primers:
        seq_str = p["sequence"]
        if seq_str not in unique_by_seq:
            unique_by_seq[seq_str] = p
        else:
            if p["match_pct"] > unique_by_seq[seq_str]["match_pct"]:
                unique_by_seq[seq_str] = p

    uniques = list(unique_by_seq.values())
    uniques.sort(key=lambda x: x["match_pct"], reverse=True)

    # 2. Filter out shorter sub-string primers when a longer parent primer is detected
    filtered_primers = []
    for p in uniques:
        is_sub = False
        for other in uniques:
            if p["primer_name"] != other["primer_name"]:
                # If p is shorter and its sequence is contained in other's sequence
                if p["sequence"] in other["sequence"] and len(p["sequence"]) < len(other["sequence"]):
                    if other["match_pct"] >= 50.0 or p["match_pct"] <= other["match_pct"] + 15.0:
                        is_sub = True
                        break
        if not is_sub:
            filtered_primers.append(p)

    return filtered_primers


def detect_primer_dimers(sequences: List[str], max_dimer_len: int = 70) -> Dict[str, Any]:
    """
    Detects potential primer dimer contamination (reads shorter than max_dimer_len bp).
    """
    if not sequences:
        return {"dimer_count": 0, "dimer_pct": 0.0, "status": "PASS"}

    n_total = len(sequences)
    dimer_count = sum(1 for s in sequences if len(s) <= max_dimer_len)
    dimer_pct = round((dimer_count / n_total) * 100, 2)

    status = "WARN" if dimer_pct > 2.0 else "PASS"
    return {
        "dimer_count": dimer_count,
        "dimer_pct": dimer_pct,
        "status": status,
        "message": f"Short fragment / primer dimer proxy (<={max_dimer_len}bp): {dimer_pct}% ({dimer_count} reads).",
    }
