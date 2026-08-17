#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd


def _first(mapping: dict[str, Any], *keys: str, default: Any = None) -> Any:
    for key in keys:
        value: Any = mapping
        for part in key.split("."):
            if not isinstance(value, dict) or part not in value:
                value = None
                break
            value = value[part]
        if value is not None:
            return value
    return default


def summarize_manifest(manifest: dict[str, Any], dataset_name: str) -> dict[str, Any]:
    """Normalize preparation/RFU run manifests into one benchmark row."""
    receptor_rows = _first(
        manifest, "receptor_row_count", "input_receptor_row_count", "input_row_count"
    )
    unique_queries = _first(
        manifest, "unique_query_count", "unique_cdr3_count", "deduplicated_query_count"
    )
    dimensions = _first(manifest, "source_atlas_dimensions", default=None)
    input_cells = _first(manifest, "input_cells", "input_cell_count", "source_cell_count")
    if input_cells is None and isinstance(dimensions, list) and dimensions:
        input_cells = dimensions[0]
    deduplication_ratio = (
        float(receptor_rows) / float(unique_queries)
        if receptor_rows is not None and unique_queries not in {None, 0}
        else None
    )
    artifact_hashes = _first(
        manifest, "rfu_artifact_hashes", "artifact_hashes", "backend_artifact_hashes", default={}
    )
    return {
        "dataset_name": dataset_name,
        "input_cells": input_cells,
        "receptor_rows": receptor_rows,
        "eligible_receptor_rows": _first(
            manifest, "eligible_receptor_row_count", "eligible_row_count"
        ),
        "unique_cdr3_queries": unique_queries,
        "deduplication_ratio": deduplication_ratio,
        "rfu_assigned_rows": _first(manifest, "assigned_row_count", "rfu_assigned_row_count"),
        "threshold_pass_rows": _first(manifest, "threshold_pass_row_count"),
        "chunk_count": _first(manifest, "chunk_count", default=0),
        "cache_reuse_count": _first(manifest, "reused_chunk_count", "cache_reuse_count", default=0),
        "elapsed_time_seconds": _first(manifest, "elapsed_seconds", "total_elapsed_seconds"),
        "peak_rss_bytes": _first(manifest, "peak_rss_bytes", "max_rss_bytes"),
        "adapter": _first(manifest, "source_adapter", "adapter", default="unknown"),
        "cache_schema_version": _first(manifest, "cache_schema_version"),
        "rfu_backend_mode": _first(manifest, "backend_mode", "mode"),
        "rfu_artifact_hashes": json.dumps(artifact_hashes, sort_keys=True),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Combine scRFU manifests into a benchmark table.")
    parser.add_argument("--manifest", type=Path, action="append", required=True)
    parser.add_argument(
        "--dataset-name",
        action="append",
        required=True,
        help="Dataset label paired positionally with --manifest (repeatable)",
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if len(args.manifest) != len(args.dataset_name):
        raise ValueError("Pass exactly one --dataset-name for each --manifest.")
    rows = []
    for path, name in zip(args.manifest, args.dataset_name, strict=True):
        rows.append(summarize_manifest(json.loads(path.read_text(encoding="utf-8")), name))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
