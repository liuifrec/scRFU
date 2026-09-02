#!/usr/bin/env python3
"""Build Figure 5 from frozen phenotype-coupling and VDJdb evidence."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from _figure_common import (
    NULL_GRAY,
    RFU_BLUE,
    SCIRPY_GREEN,
    clean_axis,
    manifest_row,
    panel_label,
    read_tsv,
    save_figure,
    set_style,
    short_cell_type,
    write_manifest,
)

SCRIPT = "manuscript/scripts/build_figure5.py"


def _phenotype_matrix(frame, top_rfus):
    selected = frame[frame["rfu_label"].isin(top_rfus)].copy()
    selected["cell_type_short"] = selected["cell_type"].map(short_cell_type)
    matrix = selected.pivot_table(
        index="cell_type_short",
        columns="rfu_label",
        values="phenotype_specific_proportion",
        aggfunc="first",
        fill_value=0,
    )
    return matrix.reindex(columns=top_rfus)


def build(args: argparse.Namespace) -> None:
    set_style()
    nearest = read_tsv(args.phenotype_nearest)
    threshold = read_tsv(args.phenotype_threshold)
    chain = read_tsv(args.native_chain_queries)
    linkage = read_tsv(args.native_chain_summary)
    sensitivity = read_tsv(args.vdjdb_sensitivity)
    nulls = read_tsv(args.vdjdb_null)

    abundance = (
        nearest[["rfu_label", "total_abundance"]]
        .drop_duplicates()
        .sort_values(["total_abundance", "rfu_label"], ascending=[False, True])
    )
    top_rfus = abundance.head(30)["rfu_label"].tolist()
    matrices = [_phenotype_matrix(nearest, top_rfus), _phenotype_matrix(threshold, top_rfus)]

    fig = plt.figure(figsize=(7.5, 6.9), constrained_layout=False)
    outer = fig.add_gridspec(
        2,
        3,
        height_ratios=[1.08, 0.92],
        left=0.16,
        right=0.95,
        top=0.91,
        bottom=0.09,
        hspace=0.42,
        wspace=0.56,
    )
    top = outer[0, :].subgridspec(1, 2, wspace=0.10)
    ax_a1 = fig.add_subplot(top[0, 0])
    ax_a2 = fig.add_subplot(top[0, 1], sharey=ax_a1)
    ax_b = fig.add_subplot(outer[1, 0])
    ax_c = fig.add_subplot(outer[1, 1])
    ax_d = fig.add_subplot(outer[1, 2])

    # A: predeclared abundance-selected phenotype profiles under both policies.
    images = []
    for ax, matrix, title in zip(
        (ax_a1, ax_a2),
        matrices,
        ("Nearest assignment", "Threshold-qualified assignment"),
        strict=True,
    ):
        images.append(ax.imshow(matrix.to_numpy(), aspect="auto", cmap="cividis", vmin=0, vmax=1))
        ax.set_title(title, loc="left", pad=4)
        ax.set_xticks(range(len(top_rfus)))
        ax.set_xticklabels(top_rfus, rotation=90, fontsize=4.4)
        ax.tick_params(length=0)
    ax_a1.set_yticks(range(len(matrices[0].index)), matrices[0].index, fontsize=5.4)
    ax_a2.tick_params(labelleft=False)
    ax_a1.set_ylabel("Cell type")
    fig.colorbar(
        images[-1], ax=ax_a2, fraction=0.045, pad=0.025, label="Within-RFU cell-type proportion"
    )
    fig.text(
        0.16,
        0.965,
        "Thirty most abundant RFUs (selection independent of phenotype)",
        fontsize=8,
        fontweight="bold",
        ha="left",
        va="top",
    )
    ax_a1.text(
        -0.16, 1.08, "A", transform=ax_a1.transAxes, fontsize=10, fontweight="bold", va="top"
    )

    # B: audited phenotype-coupling stability.
    ax_b.set_title("Phenotype-coupling stability", loc="left", pad=5)
    fractions = np.array([25, 50, 75, 100])
    cosine = np.array([0.506, 0.708, 0.863, 1.000])
    agreement = np.array([0.637, 0.735, 0.853, 1.000])
    ax_b.plot(fractions, cosine, marker="o", color=RFU_BLUE, label="Coupling cosine")
    ax_b.plot(
        fractions, agreement, marker="s", color=SCIRPY_GREEN, label="Dominant phenotype agreement"
    )
    ax_b.set_xlabel("Retained cells (%)")
    ax_b.set_ylabel("Mean stability")
    ax_b.set_xticks(fractions)
    ax_b.set_ylim(0.45, 1.02)
    ax_b.legend(frameon=False, loc="lower right")
    clean_axis(ax_b)
    panel_label(ax_b, "B")

    # C: native chain-level evidence linkage without chain expansion.
    ax_c.set_title("Native chain linkage", loc="left", pad=5)
    total_chains = len(linkage)
    eligible = int(chain["eligible"].astype(bool).sum())
    cdr3_matches = int(linkage["cdr3_has_vdjdb_evidence"].astype(bool).sum())
    cdr3_v_matches = int(linkage["cdr3_v_has_vdjdb_evidence"].astype(bool).sum())
    values = [total_chains, eligible, cdr3_matches, cdr3_v_matches]
    labels = ["AIRR chains", "Eligible TRB", "CDR3 evidence", "CDR3+V evidence"]
    colors = ["#E7EEF7", "#C8E4F2", "#F4D5A8", "#E8B66F"]
    for i, (label, value, color) in enumerate(zip(labels, values, colors, strict=True)):
        y = 0.84 - i * 0.22
        width = 0.88 - i * 0.11
        x = (1 - width) / 2
        ax_c.barh(y, width, height=0.14, left=x, color=color)
        ax_c.text(0.50, y + 0.025, label, ha="center", va="center")
        ax_c.text(0.50, y - 0.035, f"{value:,}", ha="center", va="center", fontweight="bold")
    ax_c.text(0.50, 0.02, "RFU identity: exact CDR3", ha="center", fontsize=6, color=RFU_BLUE)
    ax_c.set_xlim(0, 1)
    ax_c.set_ylim(0, 1)
    ax_c.axis("off")
    panel_label(ax_c, "C")

    # D: Wells external-annotation coherence against prespecified nulls.
    ax_d.set_title("VDJdb annotation coherence", loc="left", pad=5)
    selected_nulls = nulls[
        (nulls["dataset_label"] == "wells")
        & (nulls["match_mode"] == "cdr3")
        & (nulls["assignment_policy"] == "nearest")
        & (nulls["ambiguity_policy"] == "fractional")
        & (nulls["status"] == "completed")
    ].copy()
    names = {
        "unrestricted": "Unrestricted",
        "cdr3_length": "Length",
        "trbv": "TRBV",
        "trbv_cdr3_length": "TRBV + length",
    }
    selected_nulls["label"] = selected_nulls["null_model"].map(names)
    y = np.arange(len(selected_nulls))
    ax_d.errorbar(
        selected_nulls["null_mean"],
        y,
        xerr=selected_nulls["null_std"],
        fmt="o",
        color="#666666",
        ecolor=NULL_GRAY,
        capsize=2,
        label="Null mean ± SD",
    )
    observed = float(selected_nulls["observed"].iloc[0])
    ax_d.axvline(observed, color=RFU_BLUE, lw=1.7, label=f"Observed {observed:.3f}")
    ax_d.set_yticks(y, selected_nulls["label"])
    ax_d.invert_yaxis()
    ax_d.set_xlabel("Same-antigen pair fraction")
    ax_d.set_xlim(0.055, 0.12)
    ax_d.legend(frameon=False, loc="lower right", fontsize=5.5)
    wells_modes = sensitivity[
        (sensitivity["dataset_label"] == "wells")
        & (sensitivity["assignment_policy"] == "nearest")
        & (sensitivity["ambiguity_policy"] == "fractional")
    ].set_index("match_mode")
    ax_d.text(
        0.98,
        0.98,
        f"Exact-match coverage\nCDR3 {100 * wells_modes.loc['cdr3', 'match_fraction']:.2f}%\n"
        f"CDR3+V {100 * wells_modes.loc['cdr3_v', 'match_fraction']:.2f}%",
        transform=ax_d.transAxes,
        ha="right",
        va="top",
        fontsize=5.7,
        bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "#DDDDDD"},
    )
    ax_d.text(
        0.02,
        -0.24,
        "1,000 permutations; empirical p = 0.000999 for each null",
        transform=ax_d.transAxes,
        fontsize=5.5,
    )
    clean_axis(ax_d)
    panel_label(ax_d, "D")

    save_figure(fig, args.output_dir, "figure5")
    rows = [
        manifest_row(
            figure="Figure 5",
            panel="A",
            script=SCRIPT,
            sources=[args.phenotype_nearest, args.phenotype_threshold],
            dataset="full Wells atlas",
            analysis_unit="RFU × cell type",
            filters="top 30 RFUs by total abundance; label tie-break",
            transformation="within-RFU phenotype proportions",
            plotted_metric="cell-type proportion",
            statistical_summary="descriptive heatmap",
            comparator="nearest versus threshold assignment",
            caveat="RFUs were not selected by phenotype association; no cell-level inferential p-values.",
        ),
        manifest_row(
            figure="Figure 5",
            panel="B",
            script=SCRIPT,
            sources=[args.claim_audit],
            dataset="Wells 25k",
            analysis_unit="RFU × phenotype coupling profile",
            filters="three fixed seeds; nearest policy",
            transformation="cell subsampling relative to full bounded reference",
            plotted_metric="mean coupling cosine and dominant-phenotype agreement",
            statistical_summary="audited mean across fixed seeds",
            comparator="25%, 50%, 75%, 100% retention",
            caveat="Descriptive stability only; no cell-level inferential test.",
        ),
        manifest_row(
            figure="Figure 5",
            panel="C",
            script=SCRIPT,
            sources=[
                args.native_chain_queries,
                args.native_chain_summary,
                args.native_linkage_manifest,
            ],
            dataset="Scirpy wu2020_3k",
            analysis_unit="AIRR chain",
            filters="eligible TRB; exact CDR3 and strict chain+CDR3+normalized-V evidence keys",
            transformation="chain-aligned boolean evidence summary",
            plotted_metric="chain counts by linkage tier",
            statistical_summary="exact counts",
            comparator="CDR3 versus CDR3+V matching",
            caveat="CDR3+V evidence-query identity is distinct from CDR3-based RFU identity.",
        ),
        manifest_row(
            figure="Figure 5",
            panel="D",
            script=SCRIPT,
            sources=[args.vdjdb_sensitivity, args.vdjdb_null],
            dataset="Wells plus VDJdb 2026-06-03",
            analysis_unit="distinct matched RFU CDR3 sequence",
            filters="nearest; fractional ambiguity; CDR3 exact match",
            transformation="same-antigen pair fraction; four prespecified size-preserving nulls",
            plotted_metric="observed and null same-antigen pair fraction; match coverage",
            statistical_summary="null mean ± SD; 1,000 permutations; empirical upper-tail probability",
            comparator="unrestricted, length, TRBV, TRBV+length nulls",
            caveat="External annotation coherence is not antigen specificity or prediction.",
        ),
    ]
    write_manifest(rows, args.output_dir / "figure5_source_manifest.tsv")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--phenotype-nearest", type=Path, required=True)
    result.add_argument("--phenotype-threshold", type=Path, required=True)
    result.add_argument("--claim-audit", type=Path, required=True)
    result.add_argument("--native-chain-queries", type=Path, required=True)
    result.add_argument("--native-chain-summary", type=Path, required=True)
    result.add_argument("--native-linkage-manifest", type=Path, required=True)
    result.add_argument("--vdjdb-sensitivity", type=Path, required=True)
    result.add_argument("--vdjdb-null", type=Path, required=True)
    result.add_argument("--output-dir", type=Path, required=True)
    return result


if __name__ == "__main__":
    build(parser().parse_args())
