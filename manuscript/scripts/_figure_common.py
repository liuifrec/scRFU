"""Shared manuscript-only helpers for frozen-evidence figure prototypes."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

RFU_BLUE = "#0072B2"
THRESHOLD_BLUE = "#56B4E9"
EXACT_ORANGE = "#D89000"
SCIRPY_GREEN = "#009E73"
VJ_PURPLE = "#8E6C8A"
LENGTH_GRAY = "#7A7A7A"
DIVERSITY_BROWN = "#A67C52"
NULL_GRAY = "#B8B8B8"

REPRESENTATION_ORDER = [
    "rfu",
    "exact_cdr3",
    "scirpy_clonotype",
    "trbv_trbj",
    "cdr3_length",
    "diversity",
]
REPRESENTATION_LABELS = {
    "rfu": "RFU",
    "exact_cdr3": "Exact CDR3",
    "scirpy_clonotype": "Scirpy clonotype",
    "trbv_trbj": "TRBV+TRBJ",
    "cdr3_length": "CDR3 length",
    "diversity": "Diversity",
}
REPRESENTATION_COLORS = {
    "rfu": RFU_BLUE,
    "exact_cdr3": EXACT_ORANGE,
    "scirpy_clonotype": SCIRPY_GREEN,
    "trbv_trbj": VJ_PURPLE,
    "cdr3_length": LENGTH_GRAY,
    "diversity": DIVERSITY_BROWN,
}
DATASET_COLORS = {
    "Wells": "#4C78A8",
    "GSE190905": "#F28E2B",
    "GSE157007": "#59A14F",
    "Scirpy wu2020_3k": "#B279A2",
}

MANIFEST_COLUMNS = [
    "figure",
    "panel",
    "script",
    "source_table",
    "source_table_hash",
    "dataset",
    "analysis_unit",
    "filters",
    "transformation",
    "plotted_metric",
    "statistical_summary",
    "comparator",
    "caveat",
]


def set_style() -> None:
    """Apply the provisional manuscript style defined in figure_style.md."""
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7,
            "axes.titlesize": 8,
            "axes.labelsize": 7.5,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "axes.linewidth": 0.7,
            "lines.linewidth": 1.2,
            "lines.markersize": 4.5,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "xtick.major.size": 3,
            "ytick.major.size": 3,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def clean_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.10,
        1.06,
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
        ha="left",
    )


def sha256(path: str | Path) -> str:
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_spec(paths: list[str | Path]) -> tuple[str, str]:
    resolved = [Path(path) for path in paths]
    return ";".join(_portable_source_name(path) for path in resolved), ";".join(
        sha256(path) for path in resolved
    )


def _portable_source_name(path: Path) -> str:
    """Return an informative source label without recording a developer path."""
    parts = path.parts
    if "scRFU" in parts:
        return Path(*parts[parts.index("scRFU") + 1 :]).as_posix()
    markers = (
        "representation_consistency",
        "scirpy_comparator",
        "source_tables",
        "official_run",
        "native_vdjdb_linkage",
        "full_run",
        "downstream",
    )
    for marker in markers:
        if marker in parts:
            return Path(*parts[parts.index(marker) :]).as_posix()
    if path.name == "source_table.tsv" and len(parts) >= 2:
        return f"native_scale/{path.parent.name}/{path.name}"
    return Path(*parts[-3:]).as_posix()


def manifest_row(
    *,
    figure: str,
    panel: str,
    script: str,
    sources: list[str | Path],
    dataset: str,
    analysis_unit: str,
    filters: str,
    transformation: str,
    plotted_metric: str,
    statistical_summary: str,
    comparator: str,
    caveat: str,
) -> dict[str, Any]:
    source_names, hashes = source_spec(sources)
    return {
        "figure": figure,
        "panel": panel,
        "script": script,
        "source_table": source_names,
        "source_table_hash": hashes,
        "dataset": dataset,
        "analysis_unit": analysis_unit,
        "filters": filters,
        "transformation": transformation,
        "plotted_metric": plotted_metric,
        "statistical_summary": statistical_summary,
        "comparator": comparator,
        "caveat": caveat,
    }


def write_manifest(rows: list[dict[str, Any]], path: str | Path) -> None:
    frame = pd.DataFrame(rows, columns=MANIFEST_COLUMNS)
    frame.to_csv(path, sep="\t", index=False)


def save_figure(fig: plt.Figure, output_dir: str | Path, stem: str) -> tuple[Path, Path]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    pdf = output_dir / f"{stem}.pdf"
    png = output_dir / f"{stem}.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    return pdf, png


def read_tsv(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def short_cell_type(value: str) -> str:
    replacements = {
        "CD4-positive, alpha-beta ": "CD4 ",
        "CD8-positive, alpha-beta ": "CD8 ",
        "CD16-positive, CD56-dim natural killer cell, human": "NK (CD16+ CD56dim)",
        "CD16-negative, CD56-bright natural killer cell, human": "NK (CD56bright)",
        "T cell": "T",
        " cell": "",
        ", human": "",
    }
    output = str(value)
    for old, new in replacements.items():
        output = output.replace(old, new)
    return output
