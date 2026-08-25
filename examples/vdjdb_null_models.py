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

NULL_MODELS: dict[str, str | tuple[str, ...] | None] = {
    "unrestricted": None,
    "cdr3_length": "cdr3_length",
    "trbv": "trbv",
    "trbv_cdr3_length": ("trbv", "cdr3_length"),
}


def _read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t" if ".tsv" in path.name.lower() else ",")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_value(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, np.ndarray):
        return [_json_value(item) for item in value.tolist()]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    return value


def run_null_models(
    *,
    rfu_sequences: Path,
    evidence: Path,
    outdir: Path,
    assignment_policy: str,
    ambiguity_policy: str,
    n_permutations: int,
    random_state: int,
    save_values: bool = False,
) -> dict[str, Any]:
    sequences = _read_table(rfu_sequences)
    matches = _read_table(evidence)
    summaries: list[dict[str, Any]] = []
    values: list[pd.DataFrame] = []
    for null_model, stratification in NULL_MODELS.items():
        try:
            result = tl.rfu_antigen_permutation_test(
                sequences,
                matches,
                n_permutations=n_permutations,
                random_state=random_state,
                stratify_by=stratification,
                assignment_policy=assignment_policy,
                ambiguity_policy=ambiguity_policy,
            )
        except ValueError as exc:
            summaries.append(
                {
                    "null_model": null_model,
                    "status": "skipped",
                    "reason": str(exc),
                    "observed": np.nan,
                    "null_mean": np.nan,
                    "null_std": np.nan,
                    "empirical_upper_tail_probability": np.nan,
                    "z_score": np.nan,
                    "metric": "same_antigen_pair_fraction",
                    "n_permutations": n_permutations,
                    "random_state": random_state,
                    "stratify_by": stratification,
                    "assignment_policy": assignment_policy,
                    "ambiguity_policy": ambiguity_policy,
                }
            )
            continue
        summaries.append(
            {
                "null_model": null_model,
                "status": "completed",
                "reason": pd.NA,
                "observed": result.observed,
                "null_mean": result.null_mean,
                "null_std": result.null_std,
                "empirical_upper_tail_probability": result.empirical_upper_tail_probability,
                "z_score": result.z_score,
                **result.parameters,
            }
        )
        if save_values:
            values.append(
                pd.DataFrame(
                    {
                        "null_model": null_model,
                        "permutation_index": np.arange(n_permutations),
                        "permutation_value": result.permutation_values,
                    }
                )
            )
    outdir.mkdir(parents=True, exist_ok=True)
    summary_path = outdir / "antigen_null_model_summary.tsv"
    pd.DataFrame(summaries).to_csv(summary_path, sep="\t", index=False)
    if values:
        pd.concat(values, ignore_index=True).to_csv(
            outdir / "antigen_null_model_values.tsv.gz", sep="\t", index=False
        )
    output_files = sorted(path for path in outdir.iterdir() if path.is_file())
    manifest = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "scrfu_version": __version__,
        "parameters": {
            "assignment_policy": assignment_policy,
            "ambiguity_policy": ambiguity_policy,
            "n_permutations": n_permutations,
            "random_state": random_state,
            "null_models": list(NULL_MODELS),
            "metric": "same_antigen_pair_fraction",
            "aggregation": "unique_sequence_antigen",
        },
        "inputs": {
            "rfu_sequence_rows": len(sequences),
            "evidence_rows": len(matches),
            "rfu_sequences_sha256": _sha256(rfu_sequences),
            "evidence_sha256": _sha256(evidence),
            "runtime_only_paths": {
                "rfu_sequences": str(rfu_sequences.resolve()),
                "evidence": str(evidence.resolve()),
            },
        },
        "outputs": {
            path.name: {"sha256": _sha256(path), "size_bytes": path.stat().st_size}
            for path in output_files
        },
    }
    (outdir / "null_model_manifest.json").write_text(
        json.dumps(_json_value(manifest), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run prespecified size-preserving VDJdb antigen-coherence null models."
    )
    parser.add_argument("--rfu-sequences", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--assignment-policy", choices=("nearest", "threshold_pass"), default="nearest"
    )
    parser.add_argument(
        "--ambiguity-policy",
        choices=("fractional", "exclude_ambiguous", "multi_label"),
        default="fractional",
    )
    parser.add_argument("--n-permutations", type=int, default=1000)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--save-values", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    run_null_models(
        rfu_sequences=args.rfu_sequences,
        evidence=args.evidence,
        outdir=args.outdir,
        assignment_policy=args.assignment_policy,
        ambiguity_policy=args.ambiguity_policy,
        n_permutations=args.n_permutations,
        random_state=args.random_state,
        save_values=args.save_values,
    )


if __name__ == "__main__":
    main()
