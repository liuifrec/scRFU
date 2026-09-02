#!/usr/bin/env python3
"""Summarize RFU compression and frozen feature reuse from completed public results."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import awkward as ak
import mudata
import numpy as np
import pandas as pd

import scrfu


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _normalize_completed_rows(
    path: Path,
    *,
    dataset: str,
    sample_column: str,
    metadata_path: Path | None = None,
) -> pd.DataFrame:
    rows = pd.read_csv(path, sep="\t")
    if metadata_path is not None:
        metadata = pd.read_csv(metadata_path, sep="\t", usecols=["cell_id", sample_column])
        rows = rows.merge(metadata, on="cell_id", how="left", validate="many_to_one")
    nearest_column = "rfu_label_nearest" if "rfu_label_nearest" in rows else "rfu_label"
    pass_column = "rfu_pass_threshold" if "rfu_pass_threshold" in rows else "pass_thr"
    output = pd.DataFrame(
        {
            "dataset": dataset,
            "sample_id": rows[sample_column].astype("string"),
            "cdr3aa": rows["cdr3aa"].astype("string"),
            "rfu": rows[nearest_column].astype("string"),
            "threshold_pass": rows[pass_column].astype("boolean"),
            "rfu_score": pd.to_numeric(rows["rfu_score"], errors="coerce"),
        }
    )
    output["rfu_threshold"] = output["rfu"].where(output["threshold_pass"].fillna(False))
    output["clonotype"] = rows["clonotype_id"].astype("string") if "clonotype_id" in rows else pd.NA
    return output


def _gse190905_rows(tcr_path: Path, rfu_path: Path, clonotype_path: Path) -> pd.DataFrame:
    source = pd.read_csv(tcr_path).rename(columns={"Unnamed: 0": "cell_id"})
    source["sample_id"] = source["patient"].astype(str) + "__" + source["state"].astype(str)
    rfu = pd.read_csv(rfu_path, sep="\t")
    merged = rfu.merge(
        source[["cell_id", "sample_id"]], on="cell_id", how="left", validate="one_to_one"
    )
    output = _normalize_frame(merged, dataset="GSE190905", sample_column="sample_id")
    clonotypes = pd.read_csv(
        clonotype_path,
        sep="\t",
        usecols=["cell_id", "scirpy_clonotype"],
    ).set_index("cell_id", verify_integrity=True)
    output["clonotype"] = clonotypes.reindex(merged["cell_id"])["scirpy_clonotype"].to_numpy()
    return output


def _normalize_frame(rows: pd.DataFrame, *, dataset: str, sample_column: str) -> pd.DataFrame:
    nearest_column = "rfu_label_nearest" if "rfu_label_nearest" in rows else "rfu_label"
    pass_column = "rfu_pass_threshold" if "rfu_pass_threshold" in rows else "pass_thr"
    output = pd.DataFrame(
        {
            "dataset": dataset,
            "sample_id": rows[sample_column].astype("string"),
            "cdr3aa": rows["cdr3aa"].astype("string"),
            "rfu": rows[nearest_column].astype("string"),
            "threshold_pass": rows[pass_column].astype("boolean"),
            "rfu_score": pd.to_numeric(rows["rfu_score"], errors="coerce"),
        }
    )
    output["rfu_threshold"] = output["rfu"].where(output["threshold_pass"].fillna(False))
    output["clonotype"] = rows["clonotype_id"].astype("string") if "clonotype_id" in rows else pd.NA
    return output


def _wu_rows(path: Path) -> pd.DataFrame:
    mdata = mudata.read_h5mu(path)
    airr = mdata.mod["airr"]
    airr_records = ak.to_list(airr.obsm["airr"])
    rfu_records = ak.to_list(airr.obsm["scrfu"])
    output = []
    for cell_id, source_chains, assigned_chains in zip(
        airr.obs_names.astype(str), airr_records, rfu_records, strict=True
    ):
        if len(source_chains) != len(assigned_chains):
            raise RuntimeError("wu2020 AIRR/RFU chain alignment failed.")
        sample_id = cell_id.split("_", maxsplit=1)[0]
        clonotype = airr.obs.loc[cell_id, "clonotype_orig"]
        for source, assigned in zip(source_chains, assigned_chains, strict=True):
            if source.get("locus") != "TRB" or not assigned.get("eligible", False):
                continue
            output.append(
                {
                    "dataset": "Scirpy wu2020_3k",
                    "sample_id": sample_id,
                    "cdr3aa": source.get("junction_aa"),
                    "rfu": assigned.get("rfu_label"),
                    "threshold_pass": assigned.get("threshold_pass"),
                    "rfu_score": assigned.get("rfu_score"),
                    "clonotype": clonotype,
                }
            )
    frame = pd.DataFrame(output)
    frame["rfu_threshold"] = frame["rfu"].where(frame["threshold_pass"].fillna(False))
    return frame


def _effective_dimension(values: pd.Series) -> float:
    counts = values.dropna().value_counts().to_numpy(dtype=float)
    if not len(counts):
        return 0.0
    proportions = counts / counts.sum()
    return float(np.exp(-(proportions * np.log(proportions)).sum()))


def _feature_summary(rows: pd.DataFrame, feature: str) -> dict[str, Any]:
    valid = rows.dropna(subset=["sample_id", feature])
    sample_count = int(rows["sample_id"].nunique())
    feature_count = int(valid[feature].nunique())
    nonzero = len(valid[["sample_id", feature]].drop_duplicates())
    denominator = sample_count * feature_count
    return {
        "representation": feature,
        "features": feature_count,
        "effective_dimension": _effective_dimension(valid[feature]),
        "sample_feature_nonzero": nonzero,
        "sample_feature_sparsity": float(1 - nonzero / denominator) if denominator else np.nan,
    }


def _dataset_summary(rows: pd.DataFrame) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame]:
    sequences_per_rfu = (
        rows.dropna(subset=["cdr3aa", "rfu"])
        .drop_duplicates(["cdr3aa", "rfu"])
        .groupby("rfu", as_index=False)
        .agg(distinct_sequences=("cdr3aa", "nunique"))
    )
    prevalence = (
        rows.dropna(subset=["sample_id", "rfu"])
        .drop_duplicates(["sample_id", "rfu"])
        .groupby("rfu", as_index=False)
        .agg(samples_present=("sample_id", "nunique"))
    )
    sample_count = int(rows["sample_id"].nunique())
    prevalence["sample_prevalence"] = prevalence["samples_present"] / sample_count
    summary = {
        "dataset": str(rows["dataset"].iloc[0]),
        "receptor_rows": len(rows),
        "samples": sample_count,
        "unique_cdr3": int(rows["cdr3aa"].nunique()),
        "exact_clonotypes": (
            int(rows["clonotype"].nunique()) if rows["clonotype"].notna().any() else np.nan
        ),
        "nearest_rfus": int(rows["rfu"].nunique()),
        "threshold_qualified_rfus": int(rows["rfu_threshold"].nunique()),
        "threshold_coverage": float(rows["threshold_pass"].fillna(False).mean()),
        "median_sequences_per_rfu": float(sequences_per_rfu["distinct_sequences"].median()),
        "mean_sequences_per_rfu": float(sequences_per_rfu["distinct_sequences"].mean()),
        "maximum_sequences_per_rfu": int(sequences_per_rfu["distinct_sequences"].max()),
        "singleton_rfus": int(sequences_per_rfu["distinct_sequences"].eq(1).sum()),
        "singleton_rfu_fraction": float(sequences_per_rfu["distinct_sequences"].eq(1).mean()),
        "rfu_sample_prevalence_median": float(prevalence["sample_prevalence"].median()),
        "rfu_score_q05": float(rows["rfu_score"].quantile(0.05)),
        "rfu_score_median": float(rows["rfu_score"].median()),
        "rfu_score_q95": float(rows["rfu_score"].quantile(0.95)),
    }
    for feature in ("cdr3aa", "rfu", "rfu_threshold", "clonotype"):
        if feature == "clonotype" and not rows[feature].notna().any():
            continue
        for key, value in _feature_summary(rows, feature).items():
            if key != "representation":
                summary[f"{feature}_{key}"] = value
    return summary, sequences_per_rfu, prevalence


def _sharing(datasets: dict[str, pd.DataFrame]) -> pd.DataFrame:
    output = []
    for left_name, right_name in itertools.combinations(datasets, 2):
        for feature in ("cdr3aa", "rfu", "rfu_threshold"):
            left = set(datasets[left_name][feature].dropna().astype(str))
            right = set(datasets[right_name][feature].dropna().astype(str))
            shared = left & right
            union = left | right
            output.append(
                {
                    "dataset_left": left_name,
                    "dataset_right": right_name,
                    "representation": feature,
                    "left_features": len(left),
                    "right_features": len(right),
                    "shared_features": len(shared),
                    "jaccard": len(shared) / len(union) if union else np.nan,
                    "fraction_left_shared": len(shared) / len(left) if left else np.nan,
                    "fraction_right_shared": len(shared) / len(right) if right else np.nan,
                }
            )
    return pd.DataFrame(output)


def run(args: argparse.Namespace) -> dict[str, Any]:
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    datasets = {
        "Wells": _normalize_completed_rows(
            args.wells_results,
            dataset="Wells",
            sample_column="library_id",
            metadata_path=args.wells_metadata,
        ),
        "GSE190905": _gse190905_rows(
            args.gse190905_tcr,
            args.gse190905_results,
            args.gse190905_clonotypes,
        ),
        "GSE157007": _normalize_completed_rows(
            args.gse157007_results,
            dataset="GSE157007",
            sample_column="sample_id",
        ),
        "Scirpy wu2020_3k": _wu_rows(args.wu2020),
    }
    summaries = []
    distribution_frames = []
    prevalence_frames = []
    for name, rows in datasets.items():
        summary, distribution, prevalence = _dataset_summary(rows)
        summaries.append(summary)
        distribution.insert(0, "dataset", name)
        prevalence.insert(0, "dataset", name)
        distribution_frames.append(distribution)
        prevalence_frames.append(prevalence)
    summary_frame = pd.DataFrame(summaries)
    distribution_frame = pd.concat(distribution_frames, ignore_index=True)
    prevalence_frame = pd.concat(prevalence_frames, ignore_index=True)
    sharing_frame = _sharing(datasets)
    outputs = {
        "representation_compression.tsv": summary_frame,
        "sequences_per_rfu.tsv": distribution_frame,
        "rfu_sample_prevalence.tsv": prevalence_frame,
        "cross_dataset_feature_sharing.tsv": sharing_frame,
    }
    for filename, frame in outputs.items():
        frame.to_csv(output_dir / filename, sep="\t", index=False)
    inputs = {
        "wells_results": args.wells_results,
        "wells_metadata": args.wells_metadata,
        "gse190905_tcr": args.gse190905_tcr,
        "gse190905_results": args.gse190905_results,
        "gse190905_clonotypes": args.gse190905_clonotypes,
        "gse157007_results": args.gse157007_results,
        "wu2020": args.wu2020,
    }
    manifest = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "scrfu_version": scrfu.__version__,
        "analysis_unit": "completed receptor row; sample-feature matrices use public sample/library",
        "rfu_policy": "nearest; threshold-qualified reported separately",
        "inputs": {key: _sha256(path) for key, path in inputs.items()},
        "outputs": {filename: _sha256(output_dir / filename) for filename in outputs},
        "caveat": "Compression and feature sharing are representation properties, not biological superiority.",
    }
    (output_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--wells-results", type=Path, required=True, help="Completed Wells RFU row table."
    )
    parser.add_argument(
        "--wells-metadata",
        type=Path,
        required=True,
        help="Wells cell metadata TSV containing cell_id and library_id.",
    )
    parser.add_argument(
        "--gse190905-tcr",
        type=Path,
        required=True,
        help="Public GSE190905 receptor/metadata CSV used by the completed run.",
    )
    parser.add_argument(
        "--gse190905-results",
        type=Path,
        required=True,
        help="Completed GSE190905 RFU row table.",
    )
    parser.add_argument(
        "--gse190905-clonotypes",
        type=Path,
        required=True,
        help="Completed genuine Scirpy clonotype assignment table for GSE190905.",
    )
    parser.add_argument(
        "--gse157007-results",
        type=Path,
        required=True,
        help="Completed preregistered GSE157007 RFU row table.",
    )
    parser.add_argument(
        "--wu2020",
        type=Path,
        required=True,
        help="Public wu2020_3k H5MU with native chain-aligned scRFU annotations.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Explicit external directory for representation tables and the run manifest.",
    )
    return parser


def main() -> None:
    print(json.dumps(run(build_parser().parse_args()), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
