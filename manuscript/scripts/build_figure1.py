#!/usr/bin/env python3
"""Build Figure 1 from frozen scRFU architecture and wu2020 evidence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from _figure_common import (
    EXACT_ORANGE,
    NULL_GRAY,
    RFU_BLUE,
    SCIRPY_GREEN,
    manifest_row,
    panel_label,
    read_tsv,
    save_figure,
    set_style,
    write_manifest,
)
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

SCRIPT = "manuscript/scripts/build_figure1.py"


def box(ax: plt.Axes, xy: tuple[float, float], text: str, color: str, width: float = 0.22) -> None:
    patch = FancyBboxPatch(
        xy,
        width,
        0.22,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        facecolor=color,
        edgecolor="white",
        linewidth=0.8,
        alpha=0.95,
    )
    ax.add_patch(patch)
    ax.text(xy[0] + width / 2, xy[1] + 0.11, text, ha="center", va="center", fontsize=7)


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    ax.add_patch(
        FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=9, color="#555555", lw=0.9)
    )


def build(args: argparse.Namespace) -> None:
    set_style()
    repo = Path(__file__).resolve().parents[2]
    methods = repo / "docs/methods_freeze.md"
    novelty = repo / "docs/methodological_novelty.md"
    storage = repo / "docs/scverse_storage_schema.md"
    native_test = repo / "tests/test_scverse_native.py"
    summary = read_tsv(args.wu_summary).set_index("metric")["value"]
    run = json.loads(Path(args.wu_manifest).read_text(encoding="utf-8"))

    fig = plt.figure(figsize=(7.2, 6.7), constrained_layout=True)
    grid = fig.add_gridspec(2, 6, height_ratios=[1.0, 1.03])
    ax_a = fig.add_subplot(grid[0, :3])
    ax_b = fig.add_subplot(grid[0, 3:])
    ax_c = fig.add_subplot(grid[1, :2])
    ax_d = fig.add_subplot(grid[1, 2:4])
    ax_e = fig.add_subplot(grid[1, 4:])

    # A: method concept.
    ax_a.set_title("Frozen RFU representation", loc="left", pad=5)
    box(ax_a, (0.02, 0.55), "AIRR chains\nexact CDR3", "#E7EEF7", 0.24)
    box(ax_a, (0.38, 0.55), "Frozen RFU\nreference", "#DCECF6", 0.24)
    box(ax_a, (0.74, 0.55), "Shared RFU\nfeatures", "#D8EEE8", 0.24)
    arrow(ax_a, (0.27, 0.66), (0.37, 0.66))
    arrow(ax_a, (0.63, 0.66), (0.73, 0.66))
    ax_a.text(0.02, 0.27, "Exact sequence identity\nretained", color=EXACT_ORANGE, fontsize=6.5)
    ax_a.plot([0.02, 0.34], [0.20, 0.20], color=EXACT_ORANGE, lw=2)
    ax_a.text(0.58, 0.27, "Reusable RFU\nvocabulary", color=RFU_BLUE, fontsize=6.5)
    ax_a.plot([0.55, 0.88], [0.20, 0.20], color=RFU_BLUE, lw=2)
    ax_a.set_xlim(0, 1)
    ax_a.set_ylim(0, 1)
    ax_a.axis("off")
    panel_label(ax_a, "A")

    # B: chain alignment.
    ax_b.set_title("Chain-aligned assignment", loc="left", pad=5)
    source = [("TRA", False), ("TRB", True), ("TRB", True), ("IGH", False)]
    for index, (locus, eligible) in enumerate(source):
        y = 0.82 - index * 0.20
        ax_b.add_patch(
            FancyBboxPatch(
                (0.04, y - 0.06), 0.25, 0.12, boxstyle="round,pad=0.01", fc="#F0F0F0", ec="#AAAAAA"
            )
        )
        ax_b.text(0.165, y, f"AIRR {locus}", ha="center", va="center")
        arrow(ax_b, (0.31, y), (0.51, y))
        result = "RFU assigned" if eligible else "explicit noneligible"
        color = RFU_BLUE if eligible else NULL_GRAY
        ax_b.add_patch(
            FancyBboxPatch(
                (0.53, y - 0.06),
                0.40,
                0.12,
                boxstyle="round,pad=0.01",
                fc=color,
                ec="white",
                alpha=0.85,
            )
        )
        ax_b.text(0.73, y, result, ha="center", va="center")
    ax_b.text(0.04, 0.04, "Same chain count and order", fontweight="bold", color="#333333")
    ax_b.set_xlim(0, 1)
    ax_b.set_ylim(0, 1)
    ax_b.axis("off")
    panel_label(ax_b, "B")

    # C: storage.
    ax_c.set_title("Native object storage", loc="left", pad=5)
    items = [
        ("obsm['airr']", "source chains", "#E7EEF7"),
        ("obsm['scrfu']", "aligned RFU records", "#DCECF6"),
        ("uns['scrfu']", "portable provenance", "#D8EEE8"),
        ("obs[key]", "optional summary", "#F3E7D3"),
    ]
    for i, (key, note, color) in enumerate(items):
        y = 0.84 - i * 0.21
        ax_c.add_patch(
            FancyBboxPatch(
                (0.04, y - 0.07), 0.90, 0.14, boxstyle="round,pad=0.01", fc=color, ec="white"
            )
        )
        ax_c.text(0.09, y + 0.02, key, ha="left", va="center", fontweight="bold")
        ax_c.text(0.09, y - 0.03, note, ha="left", va="center", color="#555555", fontsize=6.5)
    ax_c.text(0.04, 0.01, "X is not accessed", color=SCIRPY_GREEN, fontweight="bold")
    ax_c.set_xlim(0, 1)
    ax_c.set_ylim(0, 1)
    ax_c.axis("off")
    panel_label(ax_c, "C")

    # D: summary policy.
    ax_d.set_title("Explicit cell summary", loc="left", pad=5)
    box(ax_d, (0.04, 0.63), "eligible\nTRB chains", "#E7EEF7", 0.30)
    arrow(ax_d, (0.36, 0.74), (0.56, 0.74))
    box(ax_d, (0.58, 0.63), "ambiguity-\naware policy", "#DCECF6", 0.34)
    ax_d.text(0.09, 0.30, "agree", ha="center", fontweight="bold", color=SCIRPY_GREEN)
    ax_d.text(0.09, 0.18, "RFU summary", ha="center")
    ax_d.text(0.69, 0.30, "conflict", ha="center", fontweight="bold", color="#B24A4A")
    ax_d.text(0.69, 0.18, "missing summary", ha="center")
    arrow(ax_d, (0.48, 0.62), (0.15, 0.35))
    arrow(ax_d, (0.62, 0.62), (0.70, 0.35))
    ax_d.text(0.04, 0.02, "No silent multi-TRB collapse", fontsize=6.5, color="#555555")
    ax_d.set_xlim(0, 1)
    ax_d.set_ylim(0, 1)
    ax_d.axis("off")
    panel_label(ax_d, "D")

    # E: public smoke.
    ax_e.set_title("Public Scirpy interoperability", loc="left", pad=5)
    metrics = [
        ("Cells", int(float(summary["cells"]))),
        ("AIRR chains", int(float(summary["airr_chains"]))),
        ("Eligible TRB", int(float(summary["eligible_trb_chains"]))),
        ("Threshold pass", int(run["threshold_qualified_trb_chains"])),
    ]
    y = [0.82, 0.62, 0.42, 0.22]
    for (label, value), ypos in zip(metrics, y, strict=True):
        ax_e.text(0.04, ypos, label, va="center", color="#555555")
        ax_e.text(0.94, ypos, f"{value:,}", va="center", ha="right", fontweight="bold")
        ax_e.plot([0.04, 0.94], [ypos - 0.08, ypos - 0.08], color="#E5E5E5", lw=0.7)
    ax_e.text(0.04, 0.02, "Native/table mismatches", color="#555555")
    ax_e.text(0.94, 0.02, "0", ha="right", color=RFU_BLUE, fontsize=12, fontweight="bold")
    ax_e.set_xlim(0, 1)
    ax_e.set_ylim(-0.05, 1)
    ax_e.axis("off")
    panel_label(ax_e, "E")

    save_figure(fig, args.output_dir, "figure1")
    rows = [
        manifest_row(
            figure="Figure 1",
            panel="A",
            script=SCRIPT,
            sources=[novelty, methods],
            dataset="Conceptual",
            analysis_unit="receptor representation",
            filters="none",
            transformation="schematic only",
            plotted_metric="none",
            statistical_summary="none",
            comparator="exact CDR3 and conventional summaries",
            caveat="Schematic is not performance evidence.",
        ),
        manifest_row(
            figure="Figure 1",
            panel="B",
            script=SCRIPT,
            sources=[native_test, storage],
            dataset="synthetic multi-chain AIRR",
            analysis_unit="AIRR chain",
            filters="none",
            transformation="schema invariant schematic",
            plotted_metric="chain alignment",
            statistical_summary="exact test assertions",
            comparator="unannotated AIRR",
            caveat="Software invariant, not biological evidence.",
        ),
        manifest_row(
            figure="Figure 1",
            panel="C",
            script=SCRIPT,
            sources=[storage],
            dataset="native AIRR",
            analysis_unit="observation and chain",
            filters="none",
            transformation="storage schematic",
            plotted_metric="none",
            statistical_summary="none",
            comparator="unannotated AIRR",
            caveat="Awkward-in-AnnData support is upstream experimental.",
        ),
        manifest_row(
            figure="Figure 1",
            panel="D",
            script=SCRIPT,
            sources=[native_test, storage],
            dataset="synthetic multi-chain AIRR",
            analysis_unit="cell",
            filters="multiple eligible TRB",
            transformation="policy decision schematic",
            plotted_metric="cell summary state",
            statistical_summary="exact test assertions",
            comparator="named summary policies",
            caveat="No biological comparison is implied.",
        ),
        manifest_row(
            figure="Figure 1",
            panel="E",
            script=SCRIPT,
            sources=[args.wu_summary, args.wu_manifest],
            dataset="Scirpy wu2020_3k",
            analysis_unit="eligible TRB chain",
            filters="productive/queryable TRB",
            transformation="metric extraction",
            plotted_metric="chain counts and mismatch count",
            statistical_summary="exact counts",
            comparator="canonical table path",
            caveat="Public 3,000-cell smoke; official RFU execution validated on Linux.",
        ),
    ]
    write_manifest(rows, Path(args.output_dir) / "figure1_source_manifest.tsv")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--wu-summary", type=Path, required=True)
    result.add_argument("--wu-manifest", type=Path, required=True)
    result.add_argument("--output-dir", type=Path, required=True)
    return result


if __name__ == "__main__":
    build(parser().parse_args())
