#!/usr/bin/env python3
"""Build Figure 3 from frozen representation-compression and sharing tables."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from _figure_common import (
    DATASET_COLORS,
    EXACT_ORANGE,
    NULL_GRAY,
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

SCRIPT = "manuscript/scripts/build_figure3.py"


def build(args: argparse.Namespace) -> None:
    set_style()
    compression = read_tsv(args.compression)
    sequences = read_tsv(args.sequences_per_rfu)
    sharing = read_tsv(args.sharing)
    comparator = read_tsv(args.representation_summary)
    datasets = compression["dataset"].tolist()

    fig = plt.figure(figsize=(7.4, 7.4), constrained_layout=False)
    grid = fig.add_gridspec(
        2,
        6,
        height_ratios=[0.96, 1.04],
        left=0.08,
        right=0.96,
        top=0.95,
        bottom=0.10,
        hspace=0.40,
        wspace=0.55,
    )
    ax_a = fig.add_subplot(grid[0, :3])
    ax_b = fig.add_subplot(grid[0, 3:])
    ax_c = fig.add_subplot(grid[1, :2])
    ax_d = fig.add_subplot(grid[1, 2:4])
    ax_e = fig.add_subplot(grid[1, 4:])

    # A: representation feature counts.
    ax_a.set_title("RFU representation compression", loc="left", pad=5)
    x = np.arange(len(datasets))
    series = [
        ("Exact CDR3", "unique_cdr3", EXACT_ORANGE),
        ("Nearest RFU", "nearest_rfus", RFU_BLUE),
        ("Threshold RFU", "threshold_qualified_rfus", THRESHOLD_BLUE),
    ]
    width = 0.23
    for offset, (label, column, color) in zip((-1, 0, 1), series, strict=True):
        ax_a.bar(x + offset * width, compression[column], width, label=label, color=color)
    ax_a.set_yscale("log")
    ax_a.set_ylabel("Observed features (log scale)")
    ax_a.set_xticks(x, [name.replace("Scirpy ", "") for name in datasets], rotation=20, ha="right")
    ax_a.legend(frameon=False, ncol=3, loc="upper right")
    clean_axis(ax_a)
    panel_label(ax_a, "A")

    # B: sequences per RFU.
    ax_b.set_title("Sequences represented per RFU", loc="left", pad=5)
    for dataset in datasets:
        values = np.sort(
            sequences.loc[sequences["dataset"] == dataset, "distinct_sequences"].to_numpy()
        )
        y = np.arange(1, len(values) + 1) / len(values)
        ax_b.step(values, y, where="post", color=DATASET_COLORS[dataset], label=dataset)
    ax_b.set_xscale("log")
    ax_b.set_xlabel("Distinct sequences per RFU (log scale)")
    ax_b.set_ylabel("Cumulative RFU fraction")
    ax_b.set_ylim(0, 1.01)
    ax_b.legend(frameon=False, loc="lower right")
    clean_axis(ax_b)
    panel_label(ax_b, "B")

    # C: feature count and sparsity, including the genuine Scirpy comparator.
    ax_c.set_title("Dimension and sparsity", loc="left", pad=5)
    display = {
        "rfu": ("RFU", RFU_BLUE),
        "exact_cdr3": ("Exact CDR3", EXACT_ORANGE),
        "scirpy_clonotype": ("Scirpy clonotype", "#009E73"),
        "trbv_trbj": ("TRBV+TRBJ", "#8E6C8A"),
        "cdr3_length": ("CDR3 length", "#7A7A7A"),
        "diversity": ("Diversity", "#A67C52"),
    }
    label_offsets = {
        "rfu": (3, 2),
        "exact_cdr3": (-46, 8),
        "scirpy_clonotype": (-50, -10),
        "trbv_trbj": (3, 2),
        "cdr3_length": (3, 2),
        "diversity": (3, -8),
    }
    for _, row in comparator.iterrows():
        label, color = display[str(row["representation"])]
        ax_c.scatter(row["features"], row["sparsity"], color=color, s=27, zorder=3)
        offset = label_offsets[str(row["representation"])]
        ax_c.annotate(
            label,
            (row["features"], row["sparsity"]),
            xytext=offset,
            textcoords="offset points",
            fontsize=5.5,
        )
    ax_c.set_xscale("log")
    ax_c.set_xlabel("Features (log scale)")
    ax_c.set_ylabel("Sample × feature sparsity")
    ax_c.set_ylim(-0.02, 1.02)
    clean_axis(ax_c)
    panel_label(ax_c, "C")

    # D: threshold coverage and vocabulary richness.
    ax_d.set_title("Frozen-reference coverage", loc="left", pad=5)
    y = np.arange(len(datasets))
    coverage = compression["threshold_coverage"].to_numpy() * 100
    ax_d.hlines(y, 70, coverage, color=NULL_GRAY, lw=2)
    ax_d.scatter(coverage, y, color=RFU_BLUE, s=30, zorder=3)
    for ypos, cov, richness in zip(
        y, coverage, compression["threshold_qualified_rfus"], strict=True
    ):
        ax_d.text(
            cov + 0.25, ypos, f"{cov:.1f}%\n{int(richness):,} RFUs", va="center", fontsize=5.5
        )
    ax_d.set_xlim(70, 82)
    ax_d.set_xlabel("Threshold-qualified sequences (%)")
    ax_d.set_yticks(y, [name.replace("Scirpy ", "") for name in datasets])
    ax_d.invert_yaxis()
    clean_axis(ax_d)
    panel_label(ax_d, "D")

    # E: exact-sequence versus RFU sharing.
    ax_e.set_title("Cross-dataset reuse", loc="left", pad=5)
    subset = sharing[sharing["representation"].isin(["cdr3aa", "rfu"])].copy()
    subset["pair"] = (
        subset["dataset_left"].str.replace("Scirpy ", "", regex=False)
        + " / "
        + subset["dataset_right"].str.replace("Scirpy ", "", regex=False)
    )
    pair_order = subset.loc[subset["representation"] == "cdr3aa", "pair"].tolist()
    matrix = subset.pivot(index="pair", columns="representation", values="jaccard").loc[
        pair_order, ["cdr3aa", "rfu"]
    ]
    image = ax_e.imshow(matrix.to_numpy(), cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax_e.set_xticks([0, 1], ["Exact CDR3", "RFU"], rotation=30, ha="right")
    ax_e.set_yticks(range(len(matrix)), matrix.index, fontsize=5.3)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix.iloc[i, j]
            ax_e.text(
                j,
                i,
                f"{value:.3f}",
                ha="center",
                va="center",
                color="white" if value > 0.55 else "#222222",
                fontsize=5.2,
            )
    fig.colorbar(image, ax=ax_e, fraction=0.05, pad=0.03, label="Jaccard")
    ax_e.tick_params(length=0)
    panel_label(ax_e, "E")

    save_figure(fig, args.output_dir, "figure3")
    rows = [
        manifest_row(
            figure="Figure 3",
            panel="A",
            script=SCRIPT,
            sources=[args.compression],
            dataset="Wells; GSE190905; GSE157007; Scirpy wu2020_3k",
            analysis_unit="dataset feature vocabulary",
            filters="canonical eligible CDR3 and frozen RFU assignments",
            transformation="log10 axis only",
            plotted_metric="unique CDR3 and RFU feature counts",
            statistical_summary="descriptive counts",
            comparator="exact CDR3; nearest RFU; threshold RFU",
            caveat="Compression is a representation property, not evidence of biological superiority.",
        ),
        manifest_row(
            figure="Figure 3",
            panel="B",
            script=SCRIPT,
            sources=[args.sequences_per_rfu],
            dataset="four public datasets",
            analysis_unit="RFU",
            filters="nearest assignment",
            transformation="within-dataset empirical cumulative distribution",
            plotted_metric="distinct CDR3 sequences per RFU",
            statistical_summary="ECDF",
            comparator="one-sequence identity baseline",
            caveat="Grouping distinct sequences does not prove functional equivalence.",
        ),
        manifest_row(
            figure="Figure 3",
            panel="C",
            script=SCRIPT,
            sources=[args.representation_summary],
            dataset="GSE190905",
            analysis_unit="sample-by-feature matrix",
            filters="identical 12 samples",
            transformation="log10 feature axis only",
            plotted_metric="feature count and sparsity",
            statistical_summary="descriptive matrix statistics",
            comparator="RFU; exact CDR3; Scirpy clonotype; TRBV+TRBJ; length; diversity",
            caveat="Lower dimensionality or sparsity is not inherently better.",
        ),
        manifest_row(
            figure="Figure 3",
            panel="D",
            script=SCRIPT,
            sources=[args.compression],
            dataset="four public datasets",
            analysis_unit="eligible receptor sequence or chain",
            filters="threshold 0.6, unchanged frozen reference",
            transformation="fraction converted to percent",
            plotted_metric="threshold coverage and qualified RFU richness",
            statistical_summary="descriptive point estimates",
            comparator="datasets under one frozen reference",
            caveat="Threshold failure is not a calibrated out-of-distribution probability.",
        ),
        manifest_row(
            figure="Figure 3",
            panel="E",
            script=SCRIPT,
            sources=[args.sharing],
            dataset="all six public dataset pairs",
            analysis_unit="dataset pair and feature identity",
            filters="nearest RFU for RFU column",
            transformation="exact set-intersection Jaccard",
            plotted_metric="cross-dataset feature Jaccard",
            statistical_summary="exact set overlap",
            comparator="exact CDR3 versus frozen RFU",
            caveat="Shared RFU feature space does not establish biological equivalence.",
        ),
    ]
    write_manifest(rows, args.output_dir / "figure3_source_manifest.tsv")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--compression", type=Path, required=True)
    result.add_argument("--sequences-per-rfu", type=Path, required=True)
    result.add_argument("--sharing", type=Path, required=True)
    result.add_argument("--representation-summary", type=Path, required=True)
    result.add_argument("--output-dir", type=Path, required=True)
    return result


if __name__ == "__main__":
    build(parser().parse_args())
