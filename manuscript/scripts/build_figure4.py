#!/usr/bin/env python3
"""Build Figure 4 from frozen GSE190905 comparator and held-out evidence."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from _figure_common import (
    REPRESENTATION_COLORS,
    REPRESENTATION_LABELS,
    REPRESENTATION_ORDER,
    RFU_BLUE,
    THRESHOLD_BLUE,
    clean_axis,
    manifest_row,
    panel_label,
    read_tsv,
    save_figure,
    set_style,
    write_manifest,
)

SCRIPT = "manuscript/scripts/build_figure4.py"


def _ordered(frame):
    return (
        frame.set_index("representation")
        .reindex(REPRESENTATION_ORDER)
        .dropna(how="all")
        .reset_index()
    )


def build(args: argparse.Namespace) -> None:
    set_style()
    pairwise = read_tsv(args.pairwise)
    retrieval = _ordered(read_tsv(args.retrieval))
    downsampling = read_tsv(args.downsampling)

    fig = plt.figure(figsize=(7.2, 6.5), constrained_layout=True)
    grid = fig.add_gridspec(2, 2)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    # A: within- versus between-donor sample similarity.
    ax_a.set_title("Repeated samples retain donor structure", loc="left", pad=5)
    cosine = pairwise[pairwise["metric"] == "cosine"]
    order = [name for name in REPRESENTATION_ORDER if name in set(cosine["representation"])]
    y = np.arange(len(order))
    for ypos, name in zip(y, order, strict=True):
        rows = cosine[cosine["representation"] == name].set_index("pair_type")
        between = float(rows.loc["between_donor", "mean"])
        within = float(rows.loc["within_donor", "mean"])
        ax_a.plot([between, within], [ypos, ypos], color="#BBBBBB", lw=1.3, zorder=1)
        ax_a.scatter(
            between,
            ypos,
            marker="o",
            facecolor="white",
            edgecolor=REPRESENTATION_COLORS[name],
            s=28,
            zorder=2,
        )
        ax_a.scatter(within, ypos, marker="o", color=REPRESENTATION_COLORS[name], s=28, zorder=3)
    ax_a.set_yticks(y, [REPRESENTATION_LABELS[name] for name in order])
    ax_a.invert_yaxis()
    ax_a.set_xlim(-0.04, 1.04)
    ax_a.set_xlabel("Mean sample-pair cosine")
    ax_a.scatter([], [], facecolor="white", edgecolor="#555555", label="Between donor")
    ax_a.scatter([], [], color="#555555", label="Within donor")
    ax_a.legend(frameon=False, loc="lower right")
    clean_axis(ax_a)
    panel_label(ax_a, "A")

    # B: donor retrieval.
    ax_b.set_title("Leave-one-timepoint-out donor retrieval", loc="left", pad=5)
    metrics = [("top1", "Top-1"), ("top3", "Top-3"), ("mean_reciprocal_rank", "MRR")]
    x = np.arange(len(retrieval))
    width = 0.23
    for offset, (column, label) in zip((-1, 0, 1), metrics, strict=True):
        ax_b.bar(
            x + offset * width,
            retrieval[column],
            width,
            label=label,
            color=["#D7E7F2", "#82B9D8", "#337FAE"][offset + 1],
        )
    ax_b.set_xticks(
        x,
        [REPRESENTATION_LABELS[name] for name in retrieval["representation"]],
        rotation=35,
        ha="right",
    )
    ax_b.set_ylim(0, 1.03)
    ax_b.set_ylabel("Retrieval score")
    ax_b.legend(frameon=False, ncol=3, loc="upper left")
    clean_axis(ax_b)
    panel_label(ax_b, "B")

    # C: downsampling stability with identical fractions and seeds.
    ax_c.set_title("Representation stability under subsampling", loc="left", pad=5)
    order = [name for name in REPRESENTATION_ORDER if name in set(downsampling["representation"])]
    for name in order:
        rows = downsampling[downsampling["representation"] == name].sort_values("fraction")
        ax_c.errorbar(
            rows["fraction"] * 100,
            rows["mean"],
            yerr=rows["std"],
            marker="o",
            capsize=2,
            color=REPRESENTATION_COLORS[name],
            label=REPRESENTATION_LABELS[name],
        )
    ax_c.set_xticks([50, 75])
    ax_c.set_xlabel("Retained receptors (%)")
    ax_c.set_ylabel("Mean cosine to full representation")
    ax_c.set_ylim(0.65, 1.01)
    ax_c.legend(frameon=False, ncol=2, loc="lower right")
    clean_axis(ax_c)
    panel_label(ax_c, "C")

    # D: preregistered held-out transfer and stability.
    ax_d.set_title("Preregistered held-out transfer (GSE157007)", loc="left", pad=5)
    card = {"transform": ax_d.transAxes, "va": "top"}
    ax_d.text(0.02, 0.98, "60,125", fontsize=11, fontweight="bold", color="#333333", **card)
    ax_d.text(0.02, 0.89, "receptors", color="#555555", **card)
    ax_d.text(0.36, 0.98, "76.88%", fontsize=11, fontweight="bold", color=RFU_BLUE, **card)
    ax_d.text(0.36, 0.89, "threshold coverage", color="#555555", fontsize=5.8, **card)
    ax_d.text(0.76, 0.98, "4,898", fontsize=11, fontweight="bold", color="#333333", **card)
    ax_d.text(0.76, 0.89, "nearest RFUs", color="#555555", **card)
    fractions = np.array([50, 75])
    nearest = np.array([0.936, 0.976])
    threshold = np.array([0.922, 0.970])
    ax_d.plot(fractions, nearest, marker="o", color=RFU_BLUE, label="Nearest RFU")
    ax_d.plot(fractions, threshold, marker="s", color=THRESHOLD_BLUE, label="Threshold RFU")
    ax_d.set_xlim(45, 80)
    ax_d.set_ylim(0.90, 0.99)
    ax_d.set_xticks(fractions)
    ax_d.set_xlabel("Retained receptors (%)")
    ax_d.set_ylabel("Mean cosine to full sample")
    ax_d.legend(frameon=False, loc="lower right")
    clean_axis(ax_d)
    panel_label(ax_d, "D")

    save_figure(fig, args.output_dir, "figure4")
    rows = [
        manifest_row(
            figure="Figure 4",
            panel="A",
            script=SCRIPT,
            sources=[args.pairwise],
            dataset="GSE190905",
            analysis_unit="patient-time sample pair",
            filters="cosine metric; same 12 samples",
            transformation="within/between donor stratification",
            plotted_metric="mean cosine similarity",
            statistical_summary="mean; source table retains SD and pair count",
            comparator="six representations in fixed order",
            caveat="Six two-visit donors; descriptive representation evidence, not deep temporal dynamics.",
        ),
        manifest_row(
            figure="Figure 4",
            panel="B",
            script=SCRIPT,
            sources=[args.retrieval],
            dataset="GSE190905",
            analysis_unit="held-out patient-time sample",
            filters="identical leave-one-timepoint-out candidates",
            transformation="none",
            plotted_metric="top-1, top-3, mean reciprocal rank",
            statistical_summary="12 queries",
            comparator="six representations in fixed order",
            caveat="RFU ties exact CDR3 on top-1 but does not lead every endpoint.",
        ),
        manifest_row(
            figure="Figure 4",
            panel="C",
            script=SCRIPT,
            sources=[args.downsampling],
            dataset="GSE190905",
            analysis_unit="sample-seed representation",
            filters="50% and 75%; fixed seeds",
            transformation="cosine to full sample vector",
            plotted_metric="mean cosine ± SD",
            statistical_summary="mean and SD across 36 observations per cell",
            comparator="six representations in fixed order",
            caveat="Coarse summaries can be more stable because they discard sequence-level detail.",
        ),
        manifest_row(
            figure="Figure 4",
            panel="D",
            script=SCRIPT,
            sources=[args.heldout_evidence, args.claim_audit, args.preregistration],
            dataset="GSE157007",
            analysis_unit="receptor and biological-sample vector",
            filters="preregistered frozen reference; three seeds",
            transformation="fraction converted to percent; audited summary extraction",
            plotted_metric="coverage, RFU richness, 50%/75% cosine stability",
            statistical_summary="descriptive counts and mean across 51 sample-seed observations per fraction",
            comparator="nearest versus threshold-qualified RFU",
            caveat="One sample per donor prevents donor retrieval; no held-out tuning was performed.",
        ),
    ]
    write_manifest(rows, args.output_dir / "figure4_source_manifest.tsv")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--pairwise", type=Path, required=True)
    result.add_argument("--retrieval", type=Path, required=True)
    result.add_argument("--downsampling", type=Path, required=True)
    result.add_argument("--heldout-evidence", type=Path, required=True)
    result.add_argument("--claim-audit", type=Path, required=True)
    result.add_argument("--preregistration", type=Path, required=True)
    result.add_argument("--output-dir", type=Path, required=True)
    return result


if __name__ == "__main__":
    build(parser().parse_args())
