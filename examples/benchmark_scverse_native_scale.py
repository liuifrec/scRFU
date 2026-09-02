#!/usr/bin/env python3
"""Benchmark bounded native AIRR assignment against a completed canonical RFU run."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import resource
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import awkward as ak
import numpy as np
import pandas as pd

import scrfu


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _text_sha256(values: list[str]) -> str:
    return hashlib.sha256("\n".join(values).encode()).hexdigest()


def _rss_kb() -> int | None:
    status = Path("/proc/self/status")
    if not status.is_file():
        return None
    for line in status.read_text(encoding="utf-8").splitlines():
        if line.startswith("VmRSS:"):
            return int(line.split()[1])
    return None


def _resource_peak_rss_kb(who: int) -> int:
    value = int(resource.getrusage(who).ru_maxrss)
    return value if platform.system() != "Darwin" else value // 1024


def _git_state() -> tuple[str | None, bool | None]:
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        dirty = bool(
            subprocess.run(
                ["git", "status", "--porcelain"],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
        )
        return commit, dirty
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None, None


def _write_timed(adata: ad.AnnData, path: Path) -> float:
    started = time.perf_counter()
    adata.write_h5ad(path)
    return time.perf_counter() - started


def _read_timed(path: Path) -> tuple[ad.AnnData, float]:
    started = time.perf_counter()
    result = ad.read_h5ad(path)
    return result, time.perf_counter() - started


def _build_airr(
    obs: pd.DataFrame,
    receptors: pd.DataFrame,
) -> ad.AnnData:
    receptor_by_cell = receptors.set_index("cell_id", verify_integrity=True)
    records: list[list[dict[str, Any]]] = []
    for cell_id in obs.index:
        if cell_id not in receptor_by_cell.index:
            records.append([])
            continue
        row = receptor_by_cell.loc[cell_id]
        records.append(
            [
                {
                    "locus": "TRB",
                    "junction_aa": row["cdr3aa"],
                    "v_call": None if pd.isna(row["v_call"]) else row["v_call"],
                    "productive": True,
                    "sequence_id": row["input_row_id"],
                }
            ]
        )
    adata = ad.AnnData(obs=obs.copy())
    adata.obsm["airr"] = ak.Array(records)
    return adata


def _parity(
    adata: ad.AnnData,
    expected: pd.DataFrame,
) -> tuple[int, float]:
    expected_by_cell = expected.set_index("cell_id", verify_integrity=True)
    records = ak.to_list(adata.obsm["scrfu"])
    mismatch_count = 0
    maximum_score_difference = 0.0
    seen = 0
    for cell_id, chains in zip(adata.obs_names.astype(str), records, strict=True):
        if not chains:
            continue
        if len(chains) != 1 or cell_id not in expected_by_cell.index:
            mismatch_count += 1
            continue
        seen += 1
        observed = chains[0]
        reference = expected_by_cell.loc[cell_id]
        observed_score = observed["rfu_score"]
        reference_score = reference["rfu_score"]
        if pd.notna(observed_score) and pd.notna(reference_score):
            maximum_score_difference = max(
                maximum_score_difference,
                abs(float(observed_score) - float(reference_score)),
            )
        equal = (
            observed["eligible"] == (reference["eligibility_status"] == "eligible")
            and observed["rfu_id"] == int(reference["rfu_id"])
            and observed["rfu_label"] == reference["rfu_label"]
            and abs(float(observed_score) - float(reference_score)) <= 1e-12
            and observed["threshold_pass"] == bool(reference["rfu_pass_threshold"])
            and observed["assignment_status"] == reference["assignment_status"]
        )
        mismatch_count += int(not equal)
    if seen != len(expected):
        mismatch_count += abs(seen - len(expected))
    return mismatch_count, maximum_score_difference


def run(args: argparse.Namespace) -> dict[str, Any]:
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    obs_path = args.completed_run / "obs_metadata.tsv.gz"
    receptors_path = args.completed_run / "receptors.tsv.gz"
    expected_path = args.completed_run / "rfu_results_per_row.tsv.gz"
    for path in (obs_path, receptors_path, expected_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    obs_source = pd.read_csv(obs_path, sep="\t", usecols=["cell_id"])
    if args.cells <= 0 or args.cells > len(obs_source):
        raise ValueError(f"--cells must be between 1 and {len(obs_source)}")
    rng = np.random.default_rng(args.seed)
    selected_positions = np.sort(rng.permutation(len(obs_source))[: args.cells])
    selected_cells = obs_source.iloc[selected_positions]["cell_id"].astype(str).tolist()
    selected = set(selected_cells)
    obs = pd.DataFrame(index=pd.Index(selected_cells, name="cell_id"))
    receptors = pd.read_csv(receptors_path, sep="\t")
    receptors = receptors.loc[receptors["cell_id"].isin(selected)].copy()
    expected = pd.read_csv(expected_path, sep="\t")
    expected = expected.loc[expected["cell_id"].isin(selected)].copy()
    if len(receptors) != len(expected):
        raise RuntimeError("Selected receptor and completed-result row counts differ.")

    rss_after_sources = _rss_kb()
    build_started = time.perf_counter()
    adata = _build_airr(obs, receptors)
    construction_seconds = time.perf_counter() - build_started
    if adata.X is not None:
        raise RuntimeError("Native scale object unexpectedly contains expression X.")
    rss_after_airr = _rss_kb()

    baseline_path = output_dir / "airr_before_scrfu.h5ad"
    baseline_write_seconds = _write_timed(adata, baseline_path)

    assignment_started = time.perf_counter()
    result = scrfu.tl.assign_rfu(
        adata,
        rfu_dir=args.rfu_dir,
        threshold=args.threshold,
        chunk_size=args.chunk_size,
        max_workers=args.max_workers,
        workdir=output_dir / "rfu_cache",
        cell_summary=False,
        inplace=False,
    )
    assignment_seconds = time.perf_counter() - assignment_started
    if result is None:
        raise RuntimeError("Native explicit assignment returned no result.")

    annotation_started = time.perf_counter()
    adata.obsm["scrfu"] = result.chain_records
    adata.uns["scrfu"] = result.provenance
    scrfu.tl.validate_scrfu_schema(adata)
    annotation_seconds = time.perf_counter() - annotation_started
    rss_after_annotation = _rss_kb()
    mismatch_count, maximum_score_difference = _parity(adata, expected)

    annotated_path = output_dir / "airr_with_chain_scrfu.h5ad"
    annotated_write_seconds = _write_timed(adata, annotated_path)
    reloaded, reload_seconds = _read_timed(annotated_path)
    reload_validation = scrfu.tl.validate_scrfu_schema(reloaded)
    reload_mismatches, reload_maximum_score_difference = _parity(reloaded, expected)

    cached_started = time.perf_counter()
    scrfu.tl.assign_rfu(
        adata,
        rfu_dir=args.rfu_dir,
        threshold=args.threshold,
        chunk_size=args.chunk_size,
        max_workers=args.max_workers,
        workdir=output_dir / "rfu_cache",
        cell_summary=True,
        summary_policy="ambiguity_aware",
        inplace=True,
    )
    cached_summary_seconds = time.perf_counter() - cached_started
    rss_after_summary = _rss_kb()
    summary_path = output_dir / "airr_with_chain_and_cell_summary.h5ad"
    summary_write_seconds = _write_timed(adata, summary_path)

    chain_counts = ak.to_numpy(ak.num(adata.obsm["airr"], axis=1))
    result_counts = ak.to_numpy(ak.num(adata.obsm["scrfu"], axis=1))
    if not np.array_equal(chain_counts, result_counts):
        raise RuntimeError("AIRR/RFU chain alignment failed.")
    if mismatch_count or reload_mismatches:
        raise RuntimeError(
            f"Native parity failed: fresh={mismatch_count}, reload={reload_mismatches}."
        )

    git_commit, git_dirty = _git_state()
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "dataset": "Wells atlas deterministic receptor-only native subset",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit,
        "git_dirty": git_dirty,
        "scrfu_version": scrfu.__version__,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "parameters": {
            "cells": args.cells,
            "seed": args.seed,
            "threshold": args.threshold,
            "chunk_size": args.chunk_size,
            "max_workers": args.max_workers,
            "sampling": "numpy default_rng permutation, source-order reconstruction",
            "chain_policy": "completed productive primary TRB; non-receptor cells retained",
        },
        "inputs": {
            "obs_metadata_sha256": _sha256(obs_path),
            "receptors_sha256": _sha256(receptors_path),
            "completed_result_sha256": _sha256(expected_path),
            "selected_cell_sha256": _text_sha256(selected_cells),
        },
        "counts": {
            "source_cells": len(obs_source),
            "selected_cells": len(obs),
            "receptor_chains": len(receptors),
            "eligible_trb": int(expected["eligibility_status"].eq("eligible").sum()),
            "unique_cdr3": int(receptors["cdr3aa"].nunique()),
            "nearest_rfus": int(expected["rfu_id"].nunique()),
            "threshold_qualified_rfus": int(
                expected.loc[expected["rfu_pass_threshold"].fillna(False), "rfu_id"].nunique()
            ),
            "threshold_qualified_chains": int(expected["rfu_pass_threshold"].fillna(False).sum()),
        },
        "threshold_coverage": float(expected["rfu_pass_threshold"].fillna(False).mean()),
        "parity": {
            "mismatch_count": mismatch_count,
            "maximum_score_absolute_difference": maximum_score_difference,
            "reload_mismatch_count": reload_mismatches,
            "reload_maximum_score_absolute_difference": reload_maximum_score_difference,
            "reload_validation": reload_validation,
        },
        "timing_seconds": {
            "airr_construction": construction_seconds,
            "baseline_serialization": baseline_write_seconds,
            "fresh_native_assignment": assignment_seconds,
            "chain_annotation_insertion": annotation_seconds,
            "annotated_serialization": annotated_write_seconds,
            "annotated_reload": reload_seconds,
            "cached_assignment_and_cell_summary": cached_summary_seconds,
            "summary_serialization": summary_write_seconds,
        },
        "memory_rss_kb": {
            "after_source_tables": rss_after_sources,
            "after_airr_construction": rss_after_airr,
            "after_chain_annotation": rss_after_annotation,
            "after_cell_summary": rss_after_summary,
            "python_parent_peak": _resource_peak_rss_kb(resource.RUSAGE_SELF),
            "rfu_child_peak": _resource_peak_rss_kb(resource.RUSAGE_CHILDREN),
        },
        "storage_bytes": {
            "airr_before_scrfu": baseline_path.stat().st_size,
            "airr_with_chain_scrfu": annotated_path.stat().st_size,
            "airr_with_chain_and_cell_summary": summary_path.stat().st_size,
            "chain_annotation_increment": annotated_path.stat().st_size
            - baseline_path.stat().st_size,
            "cell_summary_increment": summary_path.stat().st_size - annotated_path.stat().st_size,
        },
        "outputs": {
            "baseline_sha256": _sha256(baseline_path),
            "annotated_sha256": _sha256(annotated_path),
            "summary_sha256": _sha256(summary_path),
        },
        "reference_identity": result.provenance["reference_identity"],
        "artifact_hashes": result.provenance["artifact_hashes"],
        "expression_matrix": "X is None",
    }
    manifest_path = output_dir / "run_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    flat = {
        "source_cells": len(obs_source),
        "selected_cells": len(obs),
        "receptor_chains": len(receptors),
        "eligible_trb": manifest["counts"]["eligible_trb"],
        "unique_cdr3": manifest["counts"]["unique_cdr3"],
        "nearest_rfus": manifest["counts"]["nearest_rfus"],
        "threshold_qualified_rfus": manifest["counts"]["threshold_qualified_rfus"],
        "threshold_coverage": manifest["threshold_coverage"],
        "mismatch_count": mismatch_count,
        "fresh_native_assignment_seconds": assignment_seconds,
        "cached_assignment_and_summary_seconds": cached_summary_seconds,
        "python_parent_peak_rss_kb": manifest["memory_rss_kb"]["python_parent_peak"],
        "rfu_child_peak_rss_kb": manifest["memory_rss_kb"]["rfu_child_peak"],
        "baseline_size_bytes": manifest["storage_bytes"]["airr_before_scrfu"],
        "annotated_size_bytes": manifest["storage_bytes"]["airr_with_chain_scrfu"],
        "summary_size_bytes": manifest["storage_bytes"]["airr_with_chain_and_cell_summary"],
        "serialization_seconds": annotated_write_seconds,
        "reload_seconds": reload_seconds,
        "selected_cell_sha256": manifest["inputs"]["selected_cell_sha256"],
        "reference_identity": manifest["reference_identity"],
    }
    pd.DataFrame([flat]).to_csv(output_dir / "source_table.tsv", sep="\t", index=False)
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--completed-run", type=Path, required=True)
    parser.add_argument("--rfu-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cells", type=int, required=True)
    parser.add_argument("--seed", type=int, default=20260901)
    parser.add_argument("--threshold", type=float, default=0.6)
    parser.add_argument("--chunk-size", type=int, default=20000)
    parser.add_argument("--max-workers", type=int, default=2)
    return parser


def main() -> None:
    manifest = run(build_parser().parse_args())
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
