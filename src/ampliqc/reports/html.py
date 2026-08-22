"""
HTML report generator module using Jinja2 templates.
"""

from pathlib import Path
from typing import Dict, Any, Union
from jinja2 import Environment, FileSystemLoader

TEMPLATE_DIR = Path(__file__).parent / "templates"


def generate_html_report(
    sample_name: str,
    summary: Dict[str, Any],
    context_evaluation: Dict[str, Any],
    output_path: Union[str, Path]
):
    """
    Renders and saves interactive HTML report.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    env = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
    template = env.get_template("report.html.j2")

    rendered_html = template.render(
        sample_name=sample_name,
        context_name=context_evaluation.get("context", "General"),
        summary=summary,
        context_evaluation=context_evaluation,
    )

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(rendered_html)
