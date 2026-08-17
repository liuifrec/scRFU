#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib
import pandas as pd

matplotlib.use("Agg")

from scrfu import pl, tl  # noqa: E402


def _read_table(path: Path) -> pd.DataFrame:
    separator = "\t" if ".tsv" in path.name.lower() else ","
    return pd.read_csv(path, sep=separator)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _save_figure(ax: Any, path: Path) -> None:
    ax.figure.tight_layout()
    ax.figure.savefig(path, dpi=200, bbox_inches="tight")
    import matplotlib.pyplot as plt

    plt.close(ax.figure)


def run_analysis(
    *,
    rfu_results: Path,
    metadata: Path | None,
    sample_key: str,
    phenotype_key: str,
    outdir: Path,
    cell_key: str = "cell_id",
    assignment_policy: str = "nearest",
) -> dict[str, Any]:
    """Run the generic offline downstream workflow and return its manifest."""
    frame = _read_table(rfu_results)
    if cell_key not in frame:
        raise ValueError(f"RFU results are missing cell key {cell_key!r}.")
    input_rows = len(frame)
    unmatched_result_cells = 0
    unmatched_metadata_cells = 0
    if metadata is not None:
        annotations = _read_table(metadata)
        if cell_key not in annotations:
            raise ValueError(f"Metadata are missing cell key {cell_key!r}.")
        if annotations[cell_key].isna().any() or annotations[cell_key].duplicated().any():
            raise ValueError("Metadata cell identifiers must be unique and non-missing.")
        overlap = [column for column in annotations if column != cell_key and column in frame]
        if overlap:
            raise ValueError(
                "Metadata columns already present in RFU results; refusing a duplicate join: "
                f"{overlap}"
            )
        result_cells = set(frame[cell_key].dropna().astype(str))
        metadata_cells = set(annotations[cell_key].dropna().astype(str))
        unmatched_result_cells = len(result_cells - metadata_cells)
        unmatched_metadata_cells = len(metadata_cells - result_cells)
        frame = frame.merge(
            annotations, on=cell_key, how="left", sort=False, validate="many_to_one"
        )
    missing = [column for column in (sample_key, phenotype_key) if column not in frame]
    if missing:
        raise ValueError(f"Analysis input is missing requested metadata columns: {missing}")
    if frame[[sample_key, phenotype_key]].isna().any().any():
        counts = frame[[sample_key, phenotype_key]].isna().sum().to_dict()
        raise ValueError(f"Analysis metadata join left missing requested labels: {counts}")

    outdir.mkdir(parents=True, exist_ok=True)
    repertoire = tl.repertoire_metrics(frame, groupby=sample_key, weighting="cell")
    rfu_metrics = tl.rfu_metrics(
        frame,
        groupby=sample_key,
        weighting="cell",
        assignment_policy=assignment_policy,
    )
    counts = tl.rfu_pseudobulk(
        frame,
        sample_key=sample_key,
        assignment_policy=assignment_policy,
        weighting="cell",
        normalize="count",
    )
    proportions = tl.rfu_pseudobulk(
        frame,
        sample_key=sample_key,
        assignment_policy=assignment_policy,
        weighting="cell",
        normalize="proportion",
    )
    jaccard = tl.rfu_overlap(counts, metric="jaccard")
    cosine = tl.rfu_overlap(counts, metric="cosine")
    coupling = tl.rfu_phenotype_coupling(
        frame,
        phenotype_key=phenotype_key,
        sample_key=sample_key,
        assignment_policy=assignment_policy,
        weighting="cell",
    )

    tables = {
        "repertoire_metrics.tsv": repertoire,
        "rfu_metrics.tsv": rfu_metrics,
        "rfu_pseudobulk_counts.tsv": counts.matrix,
        "rfu_pseudobulk_proportions.tsv": proportions.matrix,
        "rfu_overlap_jaccard.tsv": jaccard.matrix,
        "rfu_overlap_cosine.tsv": cosine.matrix,
        "rfu_phenotype_coupling.tsv": coupling,
    }
    for name, table in tables.items():
        matrix_output = name.startswith(("rfu_pseudobulk", "rfu_overlap"))
        table.to_csv(
            outdir / name,
            sep="\t",
            index=matrix_output,
        )

    _save_figure(
        pl.rfu_metric_heatmap(proportions, value_label="Proportion"),
        outdir / "rfu_metric_heatmap.png",
    )
    _save_figure(pl.rfu_overlap_heatmap(jaccard), outdir / "rfu_overlap_heatmap.png")
    _save_figure(
        pl.rfu_convergence(rfu_metrics, annotate_top=min(5, len(rfu_metrics))),
        outdir / "rfu_convergence.png",
    )
    _save_figure(
        pl.rfu_phenotype_heatmap(coupling, phenotype_col=phenotype_key),
        outdir / "rfu_phenotype_heatmap.png",
    )

    output_paths = [outdir / name for name in tables]
    output_paths.extend(
        outdir / name
        for name in (
            "rfu_metric_heatmap.png",
            "rfu_overlap_heatmap.png",
            "rfu_convergence.png",
            "rfu_phenotype_heatmap.png",
        )
    )
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "parameters": {
            "sample_key": sample_key,
            "phenotype_key": phenotype_key,
            "cell_key": cell_key,
            "assignment_policy": assignment_policy,
            "weighting": "cell",
        },
        "inputs": {
            "rfu_results": str(rfu_results.resolve()),
            "metadata": str(metadata.resolve()) if metadata is not None else None,
            "rfu_result_rows": input_rows,
        },
        "join_qc": {
            "unmatched_result_cells": unmatched_result_cells,
            "unmatched_metadata_cells": unmatched_metadata_cells,
        },
        "outputs": {
            path.name: {"sha256": _sha256(path), "size_bytes": path.stat().st_size}
            for path in output_paths
        },
        "dimensions": {
            "pseudobulk_samples": counts.matrix.shape[0],
            "pseudobulk_rfus": counts.matrix.shape[1],
            "rfu_metric_rows": len(rfu_metrics),
            "phenotype_coupling_rows": len(coupling),
        },
    }
    (outdir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Dataset-independent scRFU downstream analysis.")
    parser.add_argument("--rfu-results", type=Path, required=True)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--sample-key", required=True)
    parser.add_argument("--phenotype-key", required=True)
    parser.add_argument("--cell-key", default="cell_id")
    parser.add_argument(
        "--assignment-policy", choices=("nearest", "threshold_pass"), default="nearest"
    )
    parser.add_argument("--outdir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    run_analysis(
        rfu_results=args.rfu_results,
        metadata=args.metadata,
        sample_key=args.sample_key,
        phenotype_key=args.phenotype_key,
        outdir=args.outdir,
        cell_key=args.cell_key,
        assignment_policy=args.assignment_policy,
    )


if __name__ == "__main__":
    main()
