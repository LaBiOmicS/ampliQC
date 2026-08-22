"""
JSON and YAML report exporter for nextQC.
"""

import json
from pathlib import Path
from typing import Dict, Any, Union
import yaml


def export_json(data: Dict[str, Any], filepath: Union[str, Path]):
    """Saves structured QC report to JSON."""
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)
    
    # Remove large non-serializable objects if present
    clean_data = dict(data)
    clean_data.pop("sequence_counts", None)

    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(clean_data, f, indent=2, ensure_ascii=False)


def export_yaml(data: Dict[str, Any], filepath: Union[str, Path]):
    """Saves structured QC report to YAML."""
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)
    
    clean_data = dict(data)
    clean_data.pop("sequence_counts", None)

    with open(filepath, "w", encoding="utf-8") as f:
        yaml.dump(clean_data, f, default_flow_style=False, sort_keys=False)
