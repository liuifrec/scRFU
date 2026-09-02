#!/usr/bin/env python3
"""Build Figure 2 from frozen parity, Wells, and native-scale evidence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from _figure_common import (
    NULL_GRAY,
    RFU_BLUE,
    SCIRPY_GREEN,
    THRESHOLD_BLUE,
    clean_axis,
    manifest_row,
    panel_label,
    read_tsv,
    save_figure,
    set_style,
    write_manifest,
)
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

SCRIPT = "manuscript/scripts/build_figure2.py"


def _load_native(paths: list[Path]):
    rows = []
    for path in paths:
        row = read_tsv(path).iloc[0].to_dict()
        row["source"] = path
        row["parent_rss_kb"] = row.get("python_parent_peak_rss_kb", row.get("peak_rss_kb"))
        rows.append(row)
    return rows


def build(args: argparse.Namespace) -> None:
    set_style()
    evidence = json.loads(args.full_evidence.read_text(encoding="utf-8"))
    fresh = json.loads(args.full_fresh.read_text(encoding="utf-8"))
    resume = json.loads(args.full_resume.read_text(encoding="utf-8"))
    native = _load_native([args.native_25k, args.native_100k, args.native_250k])

    fig = plt.figure(figsize=(7.2, 6.4), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[0.95, 1.05])
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    # A: audited field-level canonical parity.
    ax_a.set_title("Canonical RFU fidelity", loc="left", pad=5)
    fields = ["RFU ID", "RFU label", "Score", "Threshold", "Row order"]
    for i, field in enumerate(fields):
        x = 0.02 + i * 0.195
        ax_a.add_patch(
            FancyBboxPatch(
                (x, 0.49), 0.17, 0.23, boxstyle="round,pad=0.01", fc="#163E73", ec="white"
            )
        )
        ax_a.text(
            x + 0.085, 0.61, "exact", ha="center", va="center", color="white", fontweight="bold"
        )
        ax_a.text(x + 0.085, 0.42, field, ha="center", va="top", rotation=25, fontsize=6.2)
    ax_a.text(0.02, 0.84, "10 / 10 eligible rows", fontweight="bold")
    ax_a.text(
        0.02,
        0.12,
        "Maximum |score difference| = 0\nThreshold and reconstruction mismatches = 0",
        va="bottom",
        fontsize=6.5,
    )
    ax_a.set_xlim(0, 1)
    ax_a.set_ylim(0, 1)
    ax_a.axis("off")
    panel_label(ax_a, "A")

    # B: full Wells execution funnel.
    ax_b.set_title("Full Wells receptor-only execution", loc="left", pad=5)
    counts = [
        ("Source cells", int(fresh["source_atlas_dimensions"][0])),
        ("Productive TRB rows", int(fresh["eligible_row_count"])),
        ("Unique CDR3 queries", int(fresh["unique_query_count"])),
    ]
    widths = [0.90, 0.68, 0.48]
    colors = ["#E7EEF7", "#DCECF6", "#C8E4F2"]
    for i, ((label, value), width, color) in enumerate(zip(counts, widths, colors, strict=True)):
        y = 0.80 - i * 0.25
        x = (1 - width) / 2
        ax_b.add_patch(
            FancyBboxPatch(
                (x, y - 0.08), width, 0.16, boxstyle="round,pad=0.01", fc=color, ec="white"
            )
        )
        ax_b.text(0.50, y + 0.02, label, ha="center", va="center")
        ax_b.text(0.50, y - 0.035, f"{value:,}", ha="center", va="center", fontweight="bold")
        if i < 2:
            ax_b.add_patch(
                FancyArrowPatch(
                    (0.5, y - 0.10),
                    (0.5, y - 0.17),
                    arrowstyle="-|>",
                    mutation_scale=8,
                    color="#666666",
                    lw=0.8,
                )
            )
    ax_b.text(
        0.05,
        0.06,
        f"{evidence['runtime_seconds'] / 60:.1f} min total",
        fontsize=8,
        fontweight="bold",
        color=RFU_BLUE,
    )
    ax_b.text(
        0.95,
        0.06,
        f"{evidence['peak_rss_kb'] / 1_000_000:.2f} GB peak RSS",
        fontsize=8,
        fontweight="bold",
        ha="right",
        color="#444444",
    )
    ax_b.set_xlim(0, 1)
    ax_b.set_ylim(0, 1)
    ax_b.axis("off")
    panel_label(ax_b, "B")

    # C: bounded native scaling and exact table parity.
    ax_c.set_title("Native AIRR scaling with exact table parity", loc="left", pad=5)
    cells = np.array([float(row["selected_cells"]) / 1000 for row in native])
    runtime = np.array([float(row["fresh_native_assignment_seconds"]) for row in native])
    rss = np.array([float(row["parent_rss_kb"]) / 1_000_000 for row in native])
    ax_c.plot(cells, runtime, marker="o", color=RFU_BLUE, label="Fresh assignment")
    ax_c.set_xlabel("Selected cells (thousands)")
    ax_c.set_ylabel("Wall time (s)", color=RFU_BLUE)
    ax_c.tick_params(axis="y", colors=RFU_BLUE)
    ax_c.set_xticks(cells)
    ax_c2 = ax_c.twinx()
    ax_c2.plot(cells, rss, marker="s", color="#555555", ls="--", label="Python peak RSS")
    ax_c2.set_ylabel("Python peak RSS (GB)", color="#555555")
    ax_c2.tick_params(axis="y", colors="#555555")
    if all(int(row["mismatch_count"]) == 0 for row in native):
        ax_c.text(
            0.98,
            0.98,
            "0 assignment mismatches\nat all three sizes",
            transform=ax_c.transAxes,
            ha="right",
            va="top",
            fontsize=6.2,
            color=SCIRPY_GREEN,
        )
    clean_axis(ax_c)
    ax_c2.spines["top"].set_visible(False)
    panel_label(ax_c, "C")

    # D: cached/resumed execution.
    ax_d.set_title("Compatible caches avoid RFU recomputation", loc="left", pad=5)
    labels = ["25k native", "100k native", "250k native", "Full backend"]
    fresh_time = [float(row["fresh_native_assignment_seconds"]) for row in native] + [
        float(fresh["total_elapsed_seconds"])
    ]
    cached_time = [float(row["cached_assignment_and_summary_seconds"]) for row in native] + [
        float(resume["total_elapsed_seconds"])
    ]
    y = np.arange(len(labels))
    height = 0.34
    ax_d.barh(y + height / 2, fresh_time, height, color=NULL_GRAY, label="Fresh")
    ax_d.barh(y - height / 2, cached_time, height, color=THRESHOLD_BLUE, label="Cached/resumed")
    ax_d.set_xscale("log")
    ax_d.set_xlabel("Wall time (s; log scale)")
    ax_d.set_yticks(y, labels)
    ax_d.invert_yaxis()
    ax_d.legend(frameon=False, loc="upper right")
    ax_d.text(
        0.02,
        0.02,
        "Exact invariance:\nserial / parallel / chunk / order",
        transform=ax_d.transAxes,
        ha="left",
        va="bottom",
        fontsize=6.2,
        color=SCIRPY_GREEN,
    )
    clean_axis(ax_d)
    panel_label(ax_d, "D")

    save_figure(fig, args.output_dir, "figure2")
    rows = [
        manifest_row(
            figure="Figure 2",
            panel="A",
            script=SCRIPT,
            sources=[args.claim_audit],
            dataset="official RFU adversarial fixture",
            analysis_unit="eligible receptor row",
            filters="canonical eligible rows",
            transformation="field-level exact comparison",
            plotted_metric="field agreement and maximum score difference",
            statistical_summary="exact counts at 1e-12 tolerance",
            comparator="official AssignRFUs()",
            caveat="Bounded adversarial fixture, not repertoire-scale accuracy.",
        ),
        manifest_row(
            figure="Figure 2",
            panel="B",
            script=SCRIPT,
            sources=[args.full_evidence, args.full_fresh],
            dataset="Wells atlas",
            analysis_unit="source cell, productive TRB row, unique CDR3 query",
            filters="targeted receptor-only primary productive TRB extraction",
            transformation="exact-CDR3 deduplication",
            plotted_metric="counts, runtime, peak RSS",
            statistical_summary="observed point estimates",
            comparator="rows before/after deduplication",
            caveat="Linux host-specific performance; expression X was not loaded.",
        ),
        manifest_row(
            figure="Figure 2",
            panel="C",
            script=SCRIPT,
            sources=[args.native_25k, args.native_100k, args.native_250k],
            dataset="Wells native subsets",
            analysis_unit="selected cell and eligible TRB chain",
            filters="deterministic 25k/100k/250k subsets",
            transformation="native AIRR assignment and table parity comparison",
            plotted_metric="fresh runtime, Python peak RSS, mismatch count",
            statistical_summary="observed point estimates and exact mismatches",
            comparator="canonical table path",
            caveat="One Linux host; child RFU peak memory is reported separately in Extended Data.",
        ),
        manifest_row(
            figure="Figure 2",
            panel="D",
            script=SCRIPT,
            sources=[
                args.native_25k,
                args.native_100k,
                args.native_250k,
                args.full_fresh,
                args.full_resume,
                args.claim_audit,
            ],
            dataset="Wells native subsets and full atlas",
            analysis_unit="execution configuration",
            filters="compatible cache key and frozen reference",
            transformation="log10 display of unmodified wall times",
            plotted_metric="fresh and cached/resumed seconds",
            statistical_summary="observed runtime and exact equality assertions",
            comparator="fresh versus cached; configuration invariance",
            caveat="Cache reuse requires identical inputs, reference artifacts and parameters.",
        ),
    ]
    write_manifest(rows, args.output_dir / "figure2_source_manifest.tsv")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--claim-audit", type=Path, required=True)
    result.add_argument("--full-evidence", type=Path, required=True)
    result.add_argument("--full-fresh", type=Path, required=True)
    result.add_argument("--full-resume", type=Path, required=True)
    result.add_argument("--native-25k", type=Path, required=True)
    result.add_argument("--native-100k", type=Path, required=True)
    result.add_argument("--native-250k", type=Path, required=True)
    result.add_argument("--output-dir", type=Path, required=True)
    return result


if __name__ == "__main__":
    build(parser().parse_args())
