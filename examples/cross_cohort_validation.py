#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

import scrfu


def _read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t" if ".tsv" in path.name else ",")


def run(config_path: Path, outdir: Path) -> dict[str, object]:
    config = json.loads(config_path.read_text())
    reference = scrfu.tl.validate_frozen_reference(config["frozen_reference"])
    cohorts = config.get("cohorts", [])
    if not cohorts:
        raise ValueError("Configuration requires at least one cohort.")
    roles = {entry["label"]: entry["role"] for entry in cohorts}
    held_out = [label for label, role in roles.items() if role == "held_out"]
    if len(held_out) != 1:
        raise ValueError("Configuration requires exactly one held_out cohort.")
    development = [label for label, role in roles.items() if role == "development"]
    validation = [label for label, role in roles.items() if role == "validation"]
    manifest = scrfu.tl.create_heldout_validation_manifest(
        development_cohorts=development,
        validation_cohorts=validation,
        held_out_cohort=held_out[0],
        reference=reference,
        frozen_parameters=config.get("frozen_parameters", {}),
        data_hashes=config.get("data_hashes", {}),
        evaluation_metrics=config["evaluation_metrics"],
    )
    outdir.mkdir(parents=True, exist_ok=True)
    summaries: list[dict[str, object]] = []
    for entry in cohorts:
        source_path = Path(entry["results"])
        source = _read_table(source_path)
        result = scrfu.tl.transfer_cohort(
            source,
            reference,
            cohort_label=entry["label"],
            sample_key=entry["sample_key"],
            assignment_policy=entry.get("assignment_policy", "nearest"),
            weighting=entry.get("weighting", "cell"),
            observed_reference_id=entry.get("observed_reference_id"),
        )
        cohort_dir = outdir / entry["label"]
        cohort_dir.mkdir(exist_ok=True)
        result.coverage.to_csv(cohort_dir / "reference_coverage.tsv", sep="\t", index=False)
        result.score_distribution.to_csv(
            cohort_dir / "score_distribution.tsv", sep="\t", index=False
        )
        result.rfu_summary.to_csv(cohort_dir / "rfu_summary.tsv", sep="\t", index=False)
        result.sample_matrix.to_csv(cohort_dir / "rfu_pseudobulk_counts.tsv", sep="\t")
        (cohort_dir / "provenance.json").write_text(
            json.dumps(result.provenance, indent=2, sort_keys=True) + "\n"
        )
        summaries.append(
            {
                "cohort": entry["label"],
                "role": entry["role"],
                "source_rows": result.provenance["source_row_count"],
                "samples": len(result.sample_matrix),
                "rfus": result.sample_matrix.shape[1],
                "reference_verified": result.provenance["target_reference_verified"],
            }
        )
    (outdir / "heldout_manifest.json").write_text(
        json.dumps(manifest.to_dict(), indent=2, sort_keys=True) + "\n"
    )
    pd.DataFrame(summaries).to_csv(outdir / "cohort_summary.tsv", sep="\t", index=False)
    run_manifest: dict[str, object] = {
        "schema_version": "1",
        "config_path_runtime": str(config_path.resolve()),
        "frozen_reference_id": reference.immutable_reference_id,
        "cohort_count": len(summaries),
        "held_out_cohort": held_out[0],
        "no_download": True,
        "no_target_refitting": True,
    }
    (outdir / "run_manifest.json").write_text(
        json.dumps(run_manifest, indent=2, sort_keys=True) + "\n"
    )
    return run_manifest


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run configuration-frozen summaries on user-supplied local cohort results."
    )
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.outdir), sort_keys=True))


if __name__ == "__main__":
    main()
