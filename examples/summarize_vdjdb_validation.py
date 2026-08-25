#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from scrfu import __version__, tl


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _parse_dataset(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("dataset roots must use LABEL=PATH")
    label, raw_path = value.split("=", 1)
    path = Path(raw_path).expanduser().resolve()
    if not label.strip() or not path.is_dir():
        raise argparse.ArgumentTypeError(f"invalid dataset root: {value}")
    return label.strip(), path


def summarize(dataset_roots: list[tuple[str, Path]], outdir: Path) -> dict[str, Any]:
    sensitivity_rows: list[dict[str, Any]] = []
    grouping_tables: list[pd.DataFrame] = []
    null_tables: list[pd.DataFrame] = []
    score_tables: list[pd.DataFrame] = []
    input_manifests: list[dict[str, Any]] = []
    for dataset, root in dataset_roots:
        for run_dir in sorted(path for path in root.glob("*__*__*") if path.is_dir()):
            manifest_path = run_dir / "run_manifest.json"
            if not manifest_path.is_file():
                continue
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            parameters = manifest["parameters"]
            dimensions = manifest["dimensions"]
            coherence = pd.read_csv(run_dir / "rfu_antigen_coherence.tsv.gz", sep="\t")
            evidence = pd.read_csv(run_dir / "vdjdb_matches_long.tsv.gz", sep="\t")
            global_metrics = json.loads(
                (run_dir / "global_antigen_coherence.json").read_text(encoding="utf-8")
            )
            eligible = coherence.loc[coherence["eligible_for_coherence"]].copy()
            weights = eligible["coherence_contributing_sequences"].astype(float)
            weighted_normalized_entropy = (
                float(np.average(eligible["normalized_antigen_entropy"], weights=weights))
                if weights.sum()
                else np.nan
            )
            prevalence_path = run_dir / "sequence_sample_prevalence.tsv.gz"
            prevalence = (
                pd.read_csv(prevalence_path, sep="\t")["sample_prevalence_count"]
                if prevalence_path.is_file()
                else pd.Series(dtype=float)
            )
            sensitivity_rows.append(
                {
                    "dataset_label": dataset,
                    "match_mode": parameters["match_mode"],
                    "assignment_policy": parameters["assignment_policy"],
                    "ambiguity_policy": parameters["ambiguity_policy"],
                    "aggregation": "unique_sequence_antigen",
                    "biological_unit": "distinct RFU CDR3 sequence",
                    "rfu_unique_sequences": manifest["inputs"]["rfu_sequence_rows"],
                    "matched_unique_sequences": dimensions["matched_unique_sequences"],
                    "matched_query_variants": dimensions["matched_query_variants"],
                    "match_fraction": dimensions["matched_unique_sequences"]
                    / manifest["inputs"]["rfu_sequence_rows"],
                    "matched_rfus": int(coherence["vdjdb_matched_sequences"].gt(0).sum()),
                    "rfus_with_at_least_2_matched_sequences": int(
                        coherence["vdjdb_matched_sequences"].ge(2).sum()
                    ),
                    "antigen_richness": int(evidence["epitope"].dropna().nunique()),
                    "weighted_mean_rfu_antigen_purity": global_metrics[
                        "weighted_mean_rfu_antigen_purity"
                    ],
                    "weighted_mean_normalized_antigen_entropy": weighted_normalized_entropy,
                    "same_antigen_pair_fraction": global_metrics[
                        "same_antigen_pair_fraction_within_rfus"
                    ],
                    "ambiguous_matched_sequences": dimensions["ambiguous_matched_sequences"],
                    "ambiguity_fraction": dimensions["ambiguous_matched_sequence_fraction"],
                    "mean_sample_prevalence": float(prevalence.mean())
                    if len(prevalence)
                    else np.nan,
                    "median_sample_prevalence": float(prevalence.median())
                    if len(prevalence)
                    else np.nan,
                    "multi_sample_prevalence_fraction": float(prevalence.gt(1).mean())
                    if len(prevalence)
                    else np.nan,
                    "vdjdb_release": parameters["vdjdb_release"],
                    "vdjdb_sha256": manifest["reference"]["sha256"],
                    "scrfu_version": manifest["scrfu_version"],
                    "random_state": parameters["random_state"],
                    "n_permutations": parameters["n_permutations"],
                    "run_manifest_sha256": _sha256(manifest_path),
                }
            )
            sequences = pd.read_csv(manifest["inputs"]["rfu_sequences"], sep="\t")
            grouping = tl.compare_antigen_groupings(
                sequences,
                evidence,
                groupings=(
                    "rfu",
                    "trbv",
                    "cdr3_length",
                    "trbv_cdr3_length",
                    "size_matched_random",
                    "edit_distance",
                ),
                assignment_policy=parameters["assignment_policy"],
                ambiguity_policy=parameters["ambiguity_policy"],
                random_state=parameters["random_state"],
            )
            grouping.insert(0, "dataset_label", dataset)
            grouping.insert(1, "match_mode", parameters["match_mode"])
            grouping_tables.append(grouping)
            null_path = run_dir / "null_models" / "antigen_null_model_summary.tsv"
            if null_path.is_file():
                nulls = pd.read_csv(null_path, sep="\t")
                if "status" not in nulls:
                    nulls["status"] = "completed"
                    nulls["reason"] = pd.NA
                else:
                    nulls["status"] = nulls["status"].fillna("completed")
                nulls.insert(0, "dataset_label", dataset)
                nulls.insert(1, "match_mode", parameters["match_mode"])
                null_tables.append(nulls)
            scores = (
                evidence.assign(
                    evidence_score=pd.to_numeric(evidence["evidence_score"], errors="coerce")
                )
                .groupby("evidence_score", dropna=False)
                .size()
                .rename("evidence_record_count")
                .reset_index()
            )
            scores.insert(0, "dataset_label", dataset)
            scores.insert(1, "match_mode", parameters["match_mode"])
            scores.insert(2, "assignment_policy", parameters["assignment_policy"])
            scores.insert(3, "ambiguity_policy", parameters["ambiguity_policy"])
            score_tables.append(scores)
            input_manifests.append(
                {
                    "dataset_label": dataset,
                    "combination": run_dir.name,
                    "manifest_sha256": _sha256(manifest_path),
                    "reference_sha256": manifest["reference"]["sha256"],
                    "runtime_only_path": str(manifest_path),
                }
            )
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "vdjdb_sensitivity_summary.tsv": pd.DataFrame(sensitivity_rows),
        "vdjdb_grouping_comparison.tsv": pd.concat(grouping_tables, ignore_index=True),
        "vdjdb_null_summary.tsv": pd.concat(null_tables, ignore_index=True),
        "vdjdb_evidence_score_distribution.tsv": pd.concat(score_tables, ignore_index=True),
    }
    for filename, table in outputs.items():
        table.to_csv(outdir / filename, sep="\t", index=False)
    manifest = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "scrfu_version": __version__,
        "aggregation": "unique_sequence_antigen",
        "dataset_count": len(dataset_roots),
        "sensitivity_combination_count": len(sensitivity_rows),
        "inputs": input_manifests,
        "outputs": {
            filename: {
                "sha256": _sha256(outdir / filename),
                "size_bytes": (outdir / filename).stat().st_size,
                "row_count": len(table),
            }
            for filename, table in outputs.items()
        },
    }
    (outdir / "source_table_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Consolidate real VDJdb validation runs.")
    parser.add_argument("--dataset", action="append", type=_parse_dataset, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    print(json.dumps(summarize(args.dataset, args.outdir.resolve()), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
