#!/usr/bin/env python3
"""Compare frozen RFUs with genuine Scirpy clonotypes on paired public samples."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scirpy as ir

import scrfu

CHAIN_SLOTS = (
    ("TRA", "TRA_1"),
    ("TRA", "TRA_2"),
    ("TRB", "TRB_1"),
    ("TRB", "TRB_2"),
)
COUNT_REPRESENTATIONS = ("rfu", "exact_cdr3", "scirpy_clonotype", "trbv_trbj", "cdr3_length")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _chain_value(row: pd.Series, prefix: str, suffix: str) -> Any:
    value = row.get(f"{prefix}_{suffix}")
    return None if pd.isna(value) else value


def _make_scirpy_object(source: pd.DataFrame) -> Any:
    cells = []
    for _, row in source.iterrows():
        cell = ir.io.AirrCell(str(row["cell_id"]))
        for locus, prefix in CHAIN_SLOTS:
            cdr3aa = _chain_value(row, prefix, "cdr3")
            junction = _chain_value(row, prefix, "cdr3_nt")
            if cdr3aa is None and junction is None:
                continue
            chain = ir.io.AirrCell.empty_chain_dict()
            chain.update(
                {
                    "locus": locus,
                    "junction_aa": cdr3aa,
                    "junction": junction,
                    "v_call": _chain_value(row, prefix, "v_gene"),
                    "d_call": _chain_value(row, prefix, "d_gene"),
                    "j_call": _chain_value(row, prefix, "j_gene"),
                    "c_call": _chain_value(row, prefix, "c_gene"),
                    "productive": True,
                    "duplicate_count": int(_chain_value(row, prefix, "expr") or 0),
                }
            )
            cell.add_chain(chain)
        cells.append(cell)
    adata = ir.io.from_airr_cells(cells)
    for column in ("patient", "state", "method"):
        adata.obs[column] = source.set_index("cell_id").loc[adata.obs_names, column].to_numpy()
    return adata


def _count_matrix(rows: pd.DataFrame, feature: str, samples: list[str]) -> pd.DataFrame:
    valid = rows.dropna(subset=[feature]).copy()
    matrix = pd.crosstab(valid["sample_id"], valid[feature])
    return matrix.reindex(samples, fill_value=0).sort_index(axis=1)


def _cosine(left: np.ndarray, right: np.ndarray) -> float:
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    return float(np.dot(left, right) / denominator) if denominator else 0.0


def _jaccard(left: np.ndarray, right: np.ndarray) -> float:
    left_present = left > 0
    right_present = right > 0
    union = np.logical_or(left_present, right_present).sum()
    return float(np.logical_and(left_present, right_present).sum() / union) if union else 0.0


def _pairwise_summary(
    matrix: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    representation: str,
    score_name: str = "cosine",
) -> list[dict[str, Any]]:
    score = _cosine if score_name == "cosine" else _jaccard
    rows = []
    values = matrix.to_numpy(dtype=float)
    for left in range(len(matrix)):
        for right in range(left + 1, len(matrix)):
            same_donor = (
                sample_metadata.loc[matrix.index[left], "patient"]
                == sample_metadata.loc[matrix.index[right], "patient"]
            )
            rows.append(
                {
                    "representation": representation,
                    "metric": score_name,
                    "pair_type": "within_donor" if same_donor else "between_donor",
                    "value": score(values[left], values[right]),
                }
            )
    return rows


def _retrieval(
    matrix: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    representation: str,
    *,
    metric: str = "cosine",
) -> list[dict[str, Any]]:
    values = matrix.to_numpy(dtype=float)
    rows = []
    for query_index, query_id in enumerate(matrix.index):
        query_state = sample_metadata.loc[query_id, "state"]
        query_patient = sample_metadata.loc[query_id, "patient"]
        candidates = [
            index
            for index, sample_id in enumerate(matrix.index)
            if sample_metadata.loc[sample_id, "state"] != query_state
        ]
        if metric == "cosine":
            scores = [_cosine(values[query_index], values[index]) for index in candidates]
        else:
            scores = [
                1.0 / (1.0 + float(np.linalg.norm(values[query_index] - values[index])))
                for index in candidates
            ]
        ranked = sorted(
            zip(candidates, scores, strict=True),
            key=lambda item: (-item[1], matrix.index[item[0]]),
        )
        correct_rank = next(
            rank
            for rank, (index, _) in enumerate(ranked, start=1)
            if sample_metadata.loc[matrix.index[index], "patient"] == query_patient
        )
        rows.append(
            {
                "representation": representation,
                "query_sample": query_id,
                "patient": query_patient,
                "state": query_state,
                "candidate_count": len(candidates),
                "correct_rank": correct_rank,
                "top1": correct_rank <= 1,
                "top3": correct_rank <= 3,
                "reciprocal_rank": 1.0 / correct_rank,
                "similarity_metric": metric,
            }
        )
    return rows


def _diversity_matrix(rows: pd.DataFrame, samples: list[str]) -> pd.DataFrame:
    output = []
    for sample in samples:
        counts = rows.loc[rows["sample_id"].eq(sample), "exact_cdr3"].value_counts()
        proportions = counts / counts.sum()
        output.append(
            {
                "sample_id": sample,
                "richness": len(counts),
                "shannon": float(-(proportions * np.log(proportions)).sum()),
                "simpson": float(1.0 - np.square(proportions).sum()),
            }
        )
    result = pd.DataFrame(output).set_index("sample_id")
    standard_deviation = result.std(axis=0, ddof=0).replace(0, 1)
    return (result - result.mean(axis=0)) / standard_deviation


def _seed(base: int, *parts: object) -> int:
    payload = "|".join([str(base), *(str(part) for part in parts)]).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def _downsampling(
    rows: pd.DataFrame,
    matrices: dict[str, pd.DataFrame],
    *,
    seed: int,
    fractions: tuple[float, ...] = (0.5, 0.75),
    replicates: int = 3,
) -> list[dict[str, Any]]:
    output = []
    for representation in COUNT_REPRESENTATIONS:
        feature = representation
        full = matrices[representation]
        for sample_id in full.index:
            sample_rows = rows.loc[rows["sample_id"].eq(sample_id)]
            for fraction in fractions:
                size = max(1, math.floor(len(sample_rows) * fraction))
                for replicate in range(replicates):
                    rng = np.random.default_rng(
                        _seed(seed, representation, sample_id, fraction, replicate)
                    )
                    selected = sample_rows.iloc[
                        np.sort(rng.choice(len(sample_rows), size=size, replace=False))
                    ]
                    counts = selected[feature].value_counts()
                    aligned = counts.reindex(full.columns, fill_value=0).to_numpy(dtype=float)
                    output.append(
                        {
                            "representation": representation,
                            "sample_id": sample_id,
                            "fraction": fraction,
                            "replicate": replicate,
                            "sampled_cells": size,
                            "cosine_to_full": _cosine(
                                full.loc[sample_id].to_numpy(dtype=float), aligned
                            ),
                        }
                    )
    return output


def _representation_summary(
    matrices: dict[str, pd.DataFrame],
    rows: pd.DataFrame,
) -> list[dict[str, Any]]:
    output = []
    for name, matrix in matrices.items():
        values = matrix.to_numpy(dtype=float)
        if name == "diversity":
            effective_dimension = np.nan
            sparsity = np.nan
        else:
            total = values.sum()
            aggregate = values.sum(axis=0)
            proportions = aggregate[aggregate > 0] / total if total else np.array([])
            effective_dimension = (
                float(np.exp(-(proportions * np.log(proportions)).sum()))
                if len(proportions)
                else 0.0
            )
            sparsity = float(1.0 - np.count_nonzero(values) / values.size)
        output.append(
            {
                "representation": name,
                "matrix_type": (
                    "standardized sample metrics"
                    if name == "diversity"
                    else "sample-by-feature counts"
                ),
                "samples": len(matrix),
                "features": matrix.shape[1],
                "nonzero_entries": int(np.count_nonzero(values)),
                "sparsity": sparsity,
                "effective_dimension": effective_dimension,
                "represented_cells": int(rows[name].notna().sum()) if name in rows else len(matrix),
            }
        )
    return output


def run(args: argparse.Namespace) -> dict[str, Any]:
    started = time.perf_counter()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    source = pd.read_csv(args.tcr)
    source = source.rename(columns={"Unnamed: 0": "cell_id"})
    source["cell_id"] = source["cell_id"].astype(str)
    if source["cell_id"].duplicated().any():
        raise ValueError("GSE190905 source cell IDs are not unique.")

    scirpy_started = time.perf_counter()
    adata = _make_scirpy_object(source)
    ir.pp.index_chains(adata)
    ir.tl.chain_qc(adata)
    ir.pp.ir_dist(adata, metric="identity", sequence="nt", n_jobs=1)
    ir.tl.define_clonotypes(
        adata,
        key_added="scirpy_clone_id",
        receptor_arms="all",
        dual_ir="primary_only",
    )
    scirpy_seconds = time.perf_counter() - scirpy_started

    rfu = pd.read_csv(args.rfu_results, sep="\t")
    rfu = rfu.set_index("cell_id", verify_integrity=True)
    rows = source.set_index("cell_id", verify_integrity=True)
    rows["rfu"] = rfu.reindex(rows.index)["rfu_label_nearest"].astype("string")
    rows["exact_cdr3"] = rows["TRB_1_cdr3"].astype("string")
    rows["scirpy_clonotype"] = adata.obs.reindex(rows.index)["scirpy_clone_id"].astype("string")
    rows["trbv_trbj"] = (
        rows["TRB_1_v_gene"].astype("string") + "|" + rows["TRB_1_j_gene"].astype("string")
    )
    rows["cdr3_length"] = rows["exact_cdr3"].str.len().astype("Int64").astype("string")
    rows["sample_id"] = rows["patient"].astype(str) + "__" + rows["state"].astype(str)
    rows = rows.reset_index()

    sample_metadata = (
        rows[["sample_id", "patient", "state", "method"]]
        .drop_duplicates()
        .set_index("sample_id", verify_integrity=True)
        .sort_index()
    )
    samples = sample_metadata.index.tolist()
    matrices = {
        representation: _count_matrix(rows, representation, samples)
        for representation in COUNT_REPRESENTATIONS
    }
    diversity = _diversity_matrix(rows, samples)
    matrices["diversity"] = diversity

    pairwise = []
    retrieval = []
    for representation, matrix in matrices.items():
        if representation == "diversity":
            retrieval.extend(
                _retrieval(
                    matrix,
                    sample_metadata,
                    representation,
                    metric="standardized_euclidean",
                )
            )
            continue
        pairwise.extend(_pairwise_summary(matrix, sample_metadata, representation, "cosine"))
        pairwise.extend(_pairwise_summary(matrix, sample_metadata, representation, "jaccard"))
        retrieval.extend(_retrieval(matrix, sample_metadata, representation))

    pairwise_frame = pd.DataFrame(pairwise)
    pairwise_summary = (
        pairwise_frame.groupby(["representation", "metric", "pair_type"], as_index=False)["value"]
        .agg(["mean", "median", "std", "count"])
        .reset_index()
    )
    retrieval_frame = pd.DataFrame(retrieval)
    retrieval_summary = (
        retrieval_frame.groupby("representation", as_index=False)
        .agg(
            top1=("top1", "mean"),
            top3=("top3", "mean"),
            mean_reciprocal_rank=("reciprocal_rank", "mean"),
            mean_correct_rank=("correct_rank", "mean"),
            queries=("query_sample", "size"),
        )
        .sort_values("representation")
    )
    downsampling = pd.DataFrame(_downsampling(rows, matrices, seed=args.seed))
    downsampling_summary = (
        downsampling.groupby(["representation", "fraction"], as_index=False)["cosine_to_full"]
        .agg(["mean", "median", "std", "count"])
        .reset_index()
    )
    representation_summary = pd.DataFrame(_representation_summary(matrices, rows))

    assignments_path = output_dir / "scirpy_clonotype_assignments.tsv.gz"
    rows[
        [
            "cell_id",
            "patient",
            "state",
            "sample_id",
            "rfu",
            "exact_cdr3",
            "scirpy_clonotype",
            "trbv_trbj",
            "cdr3_length",
        ]
    ].to_csv(assignments_path, sep="\t", index=False, compression="gzip")
    outputs = {
        "representation_summary.tsv": representation_summary,
        "pairwise_values.tsv": pairwise_frame,
        "pairwise_summary.tsv": pairwise_summary,
        "retrieval_values.tsv": retrieval_frame,
        "retrieval_summary.tsv": retrieval_summary,
        "downsampling_values.tsv": downsampling,
        "downsampling_summary.tsv": downsampling_summary,
        "diversity_metrics_zscore.tsv": diversity.reset_index(),
    }
    for filename, frame in outputs.items():
        frame.to_csv(output_dir / filename, sep="\t", index=False)

    manifest = {
        "schema_version": 1,
        "dataset": "GSE190905",
        "role": "public paired repeated-donor representation comparator",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "scrfu_version": scrfu.__version__,
        "scirpy_version": ir.__version__,
        "python_version": platform.python_version(),
        "inputs": {
            "tcr_sha256": _sha256(args.tcr),
            "rfu_results_sha256": _sha256(args.rfu_results),
        },
        "counts": {
            "cells": len(rows),
            "patients": int(rows["patient"].nunique()),
            "samples": int(rows["sample_id"].nunique()),
            "scirpy_clonotype_assigned_cells": int(rows["scirpy_clonotype"].notna().sum()),
            "scirpy_clonotypes": int(rows["scirpy_clonotype"].nunique()),
        },
        "scirpy_parameters": {
            "index_chains_filter": ["productive", "require_junction_aa"],
            "index_chains_sort": "Scirpy 0.22.4 defaults",
            "chain_qc": True,
            "ir_dist_metric": "identity",
            "ir_dist_sequence": "nt",
            "ir_dist_n_jobs": 1,
            "define_clonotypes_receptor_arms": "all",
            "define_clonotypes_dual_ir": "primary_only",
            "define_clonotypes_partitions": "connected (function default)",
            "within_group": "receptor_type (function default)",
        },
        "comparison_parameters": {
            "biological_unit": "patient-state sample",
            "rfu_policy": "nearest frozen RFU",
            "exact_cdr3": "primary TRB_1 amino-acid CDR3",
            "vj": "primary TRB_1 V+J",
            "downsampling_fractions": [0.5, 0.75],
            "downsampling_replicates": 3,
            "random_seed": args.seed,
            "retrieval": "opposite-state candidates only; identical six-donor candidate sets",
        },
        "timing_seconds": {
            "scirpy_clonotype_construction": scirpy_seconds,
            "total": time.perf_counter() - started,
        },
        "outputs": {filename: _sha256(output_dir / filename) for filename in outputs}
        | {assignments_path.name: _sha256(assignments_path)},
        "interpretation_boundary": (
            "Scirpy clonotypes are dataset-defined exact paired receptor identities; RFUs are a "
            "frozen primary-TRB representation. Neither is expected to dominate every metric."
        ),
    }
    (output_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tcr", type=Path, required=True)
    parser.add_argument("--rfu-results", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260901)
    return parser


def main() -> None:
    print(json.dumps(run(build_parser().parse_args()), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
