#!/usr/bin/env python3
"""Build evidence-backed Extended Data prototypes 3, 5, and 6."""

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
    clean_axis,
    manifest_row,
    panel_label,
    read_tsv,
    save_figure,
    set_style,
    write_manifest,
)

SCRIPT = "manuscript/scripts/build_extended_data.py"


def _native(paths: list[Path]):
    rows = []
    for path in paths:
        row = read_tsv(path).iloc[0].to_dict()
        row["source"] = path
        row["parent_rss_kb"] = row.get("python_parent_peak_rss_kb", row.get("peak_rss_kb"))
        row["child_rss_kb"] = row.get("rfu_child_peak_rss_kb", np.nan)
        rows.append(row)
    return rows


def build_ed3(args: argparse.Namespace, rows: list[dict]) -> list[dict]:
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0), constrained_layout=True)
    cells = np.array([float(row["selected_cells"]) / 1000 for row in rows])
    labels = [f"{int(value)}k" for value in cells]

    ax = axes[0, 0]
    states = [
        ("AIRR baseline", "baseline_size_bytes", "#D9E4F2"),
        ("Chain annotation", "annotated_size_bytes", "#86B8D7"),
        ("Cell summary", "summary_size_bytes", RFU_BLUE),
    ]
    x = np.arange(len(rows))
    width = 0.24
    for offset, (label, field, color) in zip((-1, 0, 1), states, strict=True):
        ax.bar(
            x + offset * width,
            [row[field] / 1_000_000 for row in rows],
            width,
            label=label,
            color=color,
        )
    ax.set_xticks(x, labels)
    ax.set_xlabel("Selected cells")
    ax.set_ylabel("Serialized object size (MB)")
    ax.legend(frameon=False)
    ax.set_title("Native storage overhead", loc="left")
    clean_axis(ax)
    panel_label(ax, "A")

    ax = axes[0, 1]
    ax.plot(
        cells,
        [row["serialization_seconds"] for row in rows],
        marker="o",
        color=RFU_BLUE,
        label="Write",
    )
    ax.plot(
        cells, [row["reload_seconds"] for row in rows], marker="s", color="#666666", label="Reload"
    )
    ax.set_xlabel("Selected cells (thousands)")
    ax.set_ylabel("Wall time (s)")
    ax.set_title("H5AD write and reload", loc="left")
    ax.legend(frameon=False)
    clean_axis(ax)
    panel_label(ax, "B")

    ax = axes[1, 0]
    ax.plot(
        cells,
        [row["parent_rss_kb"] / 1_000_000 for row in rows],
        marker="o",
        color=RFU_BLUE,
        label="Python parent",
    )
    child = np.array([row["child_rss_kb"] / 1_000_000 for row in rows], dtype=float)
    mask = np.isfinite(child)
    ax.plot(cells[mask], child[mask], marker="s", color="#666666", ls="--", label="RFU child")
    ax.set_xlabel("Selected cells (thousands)")
    ax.set_ylabel("Peak RSS (GB)")
    ax.set_title("Process memory", loc="left")
    ax.legend(frameon=False)
    clean_axis(ax)
    panel_label(ax, "C")

    ax = axes[1, 1]
    checks = [
        "H5AD round trip",
        "H5MU round trip",
        "Subset",
        "Compatible concat",
        "Incompatible reject",
        "X=None",
    ]
    matrix = np.ones((1, len(checks)))
    ax.imshow(matrix, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(checks)), checks, rotation=35, ha="right")
    ax.set_yticks([0], ["Exact invariant"])
    for i in range(len(checks)):
        ax.text(
            i, 0, "PASS", color="white", ha="center", va="center", fontweight="bold", fontsize=6
        )
    ax.tick_params(length=0)
    ax.set_title("Standard scverse operations", loc="left")
    panel_label(ax, "D")

    save_figure(fig, args.output_dir, "extended_data3")
    common = {
        "figure": "Extended Data 3",
        "script": SCRIPT,
        "sources": [row["source"] for row in rows],
        "dataset": "Wells native 25k/100k/250k",
        "analysis_unit": "receptor-only native object",
        "filters": "deterministic subsets",
        "statistical_summary": "observed point estimates",
        "caveat": "Storage and memory depend on host and HDF5 compression.",
    }
    return [
        manifest_row(
            panel="A",
            transformation="bytes to decimal MB",
            plotted_metric="serialized object size",
            comparator="AIRR baseline; chain annotation; cell summary",
            **common,
        ),
        manifest_row(
            panel="B",
            transformation="none",
            plotted_metric="write and reload wall time",
            comparator="write versus reload",
            **common,
        ),
        manifest_row(
            panel="C",
            transformation="kB to decimal GB",
            plotted_metric="peak RSS",
            comparator="Python parent versus RFU child process",
            **common,
        ),
        manifest_row(
            figure="Extended Data 3",
            panel="D",
            script=SCRIPT,
            sources=[args.native_test, args.storage_schema],
            dataset="synthetic native fixtures",
            analysis_unit="observation and AIRR chain",
            filters="supported schema",
            transformation="technical assertion matrix",
            plotted_metric="pass/reject status",
            statistical_summary="exact tests",
            comparator="before and after scverse operation",
            caveat="Technical interoperability evidence, not a biological effect.",
        ),
    ]


def build_ed5(args: argparse.Namespace) -> list[dict]:
    representation = read_tsv(args.representation_summary).set_index("representation")
    retrieval = read_tsv(args.retrieval).set_index("representation")
    pairwise = read_tsv(args.pairwise)
    downsampling = read_tsv(args.downsampling)
    order = [name for name in REPRESENTATION_ORDER if name in representation.index]
    labels = [REPRESENTATION_LABELS[name] for name in order]
    colors = [REPRESENTATION_COLORS[name] for name in order]

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0), constrained_layout=True)
    ax = axes[0, 0]
    ax.scatter(
        representation.loc[order, "features"], representation.loc[order, "sparsity"], c=colors, s=35
    )
    label_offsets = {
        "rfu": (3, 2),
        "exact_cdr3": (-48, 8),
        "scirpy_clonotype": (3, -10),
        "trbv_trbj": (3, 2),
        "cdr3_length": (3, 2),
        "diversity": (3, -8),
    }
    for name in order:
        ax.annotate(
            REPRESENTATION_LABELS[name],
            (representation.loc[name, "features"], representation.loc[name, "sparsity"]),
            xytext=label_offsets[name],
            textcoords="offset points",
            fontsize=5.5,
        )
    ax.set_xscale("log")
    ax.set_xlabel("Features (log scale)")
    ax.set_ylabel("Sparsity")
    ax.set_title("Representation scale", loc="left")
    clean_axis(ax)
    panel_label(ax, "A")

    ax = axes[0, 1]
    cosine = pairwise[pairwise["metric"] == "cosine"].pivot(
        index="representation", columns="pair_type", values="mean"
    )
    pair_order = [name for name in order if name in cosine.index]
    difference = cosine.loc[pair_order, "within_donor"] - cosine.loc[pair_order, "between_donor"]
    ax.barh(
        np.arange(len(pair_order)),
        difference,
        color=[REPRESENTATION_COLORS[name] for name in pair_order],
    )
    ax.set_yticks(range(len(pair_order)), [REPRESENTATION_LABELS[name] for name in pair_order])
    ax.invert_yaxis()
    ax.axvline(0, color="#333333", lw=0.7)
    ax.set_xlabel("Within − between donor mean cosine")
    ax.set_title("Donor contrast", loc="left")
    clean_axis(ax)
    panel_label(ax, "B")

    ax = axes[1, 0]
    width = 0.36
    x = np.arange(len(order))
    ax.bar(x - width / 2, retrieval.loc[order, "top1"], width, color="#9ECAE1", label="Top-1")
    ax.bar(
        x + width / 2,
        retrieval.loc[order, "mean_reciprocal_rank"],
        width,
        color="#3182BD",
        label="MRR",
    )
    ax.set_xticks(x, labels, rotation=35, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("Retrieval score")
    ax.set_title("Donor retrieval", loc="left")
    ax.legend(frameon=False)
    clean_axis(ax)
    panel_label(ax, "C")

    ax = axes[1, 1]
    for name in order:
        rows = downsampling[downsampling["representation"] == name].sort_values("fraction")
        ax.errorbar(
            rows["fraction"] * 100,
            rows["mean"],
            yerr=rows["std"],
            marker="o",
            capsize=2,
            color=REPRESENTATION_COLORS[name],
            label=REPRESENTATION_LABELS[name],
        )
    ax.set_xlabel("Retained receptors (%)")
    ax.set_ylabel("Cosine to full representation")
    ax.set_title("Subsampling stability", loc="left")
    ax.legend(frameon=False, ncol=2, fontsize=5.5)
    clean_axis(ax)
    panel_label(ax, "D")

    save_figure(fig, args.output_dir, "extended_data5")
    common = {
        "figure": "Extended Data 5",
        "script": SCRIPT,
        "dataset": "GSE190905",
        "filters": "identical 12 samples and candidate sets",
        "caveat": "Different representations lead different endpoints; RFU is not universally superior.",
    }
    return [
        manifest_row(
            panel="A",
            sources=[args.representation_summary],
            analysis_unit="sample-by-feature matrix",
            transformation="log10 feature axis only",
            plotted_metric="feature count and sparsity",
            statistical_summary="descriptive matrix statistics",
            comparator="six representations",
            **common,
        ),
        manifest_row(
            panel="B",
            sources=[args.pairwise],
            analysis_unit="patient-time sample pair",
            transformation="within mean minus between mean",
            plotted_metric="cosine donor contrast",
            statistical_summary="difference of source means",
            comparator="five defined representations",
            **common,
        ),
        manifest_row(
            panel="C",
            sources=[args.retrieval],
            analysis_unit="held-out sample",
            transformation="none",
            plotted_metric="top-1 and MRR",
            statistical_summary="12 queries",
            comparator="six representations",
            **common,
        ),
        manifest_row(
            panel="D",
            sources=[args.downsampling],
            analysis_unit="sample-seed representation",
            transformation="cosine to full vector",
            plotted_metric="mean cosine ± SD",
            statistical_summary="36 observations per cell",
            comparator="six representations",
            **common,
        ),
    ]


def build_ed6(args: argparse.Namespace) -> list[dict]:
    sensitivity = read_tsv(args.vdjdb_sensitivity)
    grouping = read_tsv(args.vdjdb_grouping)
    nulls = read_tsv(args.vdjdb_null)
    datasets = ["wells", "GSE190905", "GSE157007"]
    configs = sensitivity[["match_mode", "assignment_policy", "ambiguity_policy"]].drop_duplicates()
    configs["label"] = (
        configs["match_mode"]
        + " | "
        + configs["assignment_policy"]
        + " | "
        + configs["ambiguity_policy"]
    )

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.6), constrained_layout=True)
    for ax, metric, title, fmt in [
        (axes[0, 0], "match_fraction", "Exact-match coverage", ".3f"),
        (axes[0, 1], "ambiguity_fraction", "Annotation ambiguity", ".2f"),
    ]:
        joined = configs.merge(
            sensitivity, on=["match_mode", "assignment_policy", "ambiguity_policy"]
        )
        matrix = joined.pivot(index="label", columns="dataset_label", values=metric).reindex(
            index=configs["label"], columns=datasets
        )
        image = ax.imshow(
            matrix.to_numpy(),
            cmap="Blues",
            vmin=0,
            vmax=np.nanmax(matrix.to_numpy()),
            aspect="auto",
        )
        ax.set_xticks(
            range(len(datasets)), ["Wells", "GSE190905", "GSE157007"], rotation=25, ha="right"
        )
        ax.set_yticks(
            range(len(matrix)),
            [
                value.replace("threshold_pass", "threshold").replace("exclude_ambiguous", "exclude")
                for value in matrix.index
            ],
            fontsize=5.1,
        )
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                value = matrix.iloc[i, j]
                ax.text(
                    j,
                    i,
                    format(value, fmt),
                    ha="center",
                    va="center",
                    fontsize=5.2,
                    color="white" if value > 0.65 * np.nanmax(matrix.to_numpy()) else "#222222",
                )
        fig.colorbar(image, ax=ax, fraction=0.04, pad=0.02)
        ax.set_title(title, loc="left")
        ax.tick_params(length=0)
    panel_label(axes[0, 0], "A")
    panel_label(axes[0, 1], "B")

    ax = axes[1, 0]
    selected = grouping[
        (grouping["dataset_label"] == "wells")
        & (grouping["match_mode"] == "cdr3")
        & (grouping["assignment_policy"] == "nearest")
        & (grouping["ambiguity_policy"] == "fractional")
        & (grouping["metric"] == "same_antigen_pair_fraction")
        & (grouping["status"] == "completed")
    ].copy()
    methods = [
        "rfu",
        "trbv",
        "cdr3_length",
        "trbv_cdr3_length",
        "edit_distance",
        "size_matched_random",
    ]
    selected = (
        selected.set_index("grouping_method")
        .reindex(methods)
        .dropna(subset=["value"])
        .reset_index()
    )
    colors = [RFU_BLUE] + ["#777777"] * (len(selected) - 1)
    ax.barh(np.arange(len(selected)), selected["value"], color=colors)
    ax.set_yticks(range(len(selected)), selected["grouping_method"].str.replace("_", " "))
    ax.invert_yaxis()
    ax.set_xlabel("Same-antigen pair fraction")
    ax.set_title("Grouping baselines (Wells, CDR3)", loc="left")
    clean_axis(ax)
    panel_label(ax, "C")

    ax = axes[1, 1]
    selected = nulls[
        (nulls["match_mode"] == "cdr3")
        & (nulls["assignment_policy"] == "nearest")
        & (nulls["ambiguity_policy"] == "fractional")
        & (nulls["status"] == "completed")
    ].copy()
    for dataset in datasets:
        rows = selected[selected["dataset_label"] == dataset]
        ax.scatter(
            rows["null_model"]
            .str.replace("trbv_cdr3_length", "TRBV+length")
            .str.replace("cdr3_length", "Length")
            .str.replace("trbv", "TRBV")
            .str.replace("unrestricted", "Unrestricted"),
            rows["observed"] - rows["null_mean"],
            label=dataset.replace("wells", "Wells"),
            s=24,
        )
    ax.axhline(0, color="#555555", lw=0.7)
    ax.set_ylabel("Observed − null pair fraction")
    ax.set_title("Prespecified null sensitivity", loc="left")
    ax.tick_params(axis="x", rotation=30)
    ax.legend(frameon=False)
    clean_axis(ax)
    panel_label(ax, "D")

    save_figure(fig, args.output_dir, "extended_data6")
    sensitivity_common = {
        "figure": "Extended Data 6",
        "script": SCRIPT,
        "sources": [args.vdjdb_sensitivity],
        "dataset": "Wells; GSE190905; GSE157007 plus VDJdb",
        "analysis_unit": "distinct matched sequence",
        "filters": "all 24 frozen match/assignment/ambiguity combinations",
        "transformation": "policy matrix",
        "statistical_summary": "descriptive frozen combinations",
        "comparator": "matching and ambiguity policies",
        "caveat": "Sparse exact matches remain visible and were not optimized away.",
    }
    return [
        manifest_row(panel="A", plotted_metric="exact-match fraction", **sensitivity_common),
        manifest_row(panel="B", plotted_metric="ambiguity fraction", **sensitivity_common),
        manifest_row(
            figure="Extended Data 6",
            panel="C",
            script=SCRIPT,
            sources=[args.vdjdb_grouping],
            dataset="Wells plus VDJdb",
            analysis_unit="sequence grouping",
            filters="CDR3; nearest; fractional ambiguity",
            transformation="none",
            plotted_metric="same-antigen pair fraction",
            statistical_summary="descriptive grouping comparison",
            comparator="RFU and five receptor-property/random groupings",
            caveat="External annotation coherence is not antigen specificity.",
        ),
        manifest_row(
            figure="Extended Data 6",
            panel="D",
            script=SCRIPT,
            sources=[args.vdjdb_null],
            dataset="three public datasets plus VDJdb",
            analysis_unit="matched label assignment",
            filters="CDR3; nearest; fractional; completed 1,000-permutation cells",
            transformation="observed minus null mean",
            plotted_metric="same-antigen pair-fraction difference",
            statistical_summary="four size-preserving nulls",
            comparator="unrestricted, length, TRBV, TRBV+length",
            caveat="Nulls control only named receptor properties.",
        ),
    ]


def build(args: argparse.Namespace) -> None:
    set_style()
    rows = _native([args.native_25k, args.native_100k, args.native_250k])
    manifest = build_ed3(args, rows) + build_ed5(args) + build_ed6(args)
    write_manifest(manifest, args.output_dir / "extended_data_source_manifest.tsv")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    for name in (
        "native-25k",
        "native-100k",
        "native-250k",
        "native-test",
        "storage-schema",
        "representation-summary",
        "retrieval",
        "pairwise",
        "downsampling",
        "vdjdb-sensitivity",
        "vdjdb-grouping",
        "vdjdb-null",
    ):
        result.add_argument(f"--{name}", type=Path, required=True)
    result.add_argument("--output-dir", type=Path, required=True)
    return result


if __name__ == "__main__":
    build(parser().parse_args())
