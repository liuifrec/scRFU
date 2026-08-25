"""Validation of completed table-level RFU runs without backend recomputation."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .io import file_sha256

_REQUIRED_FILES = (
    "receptors.tsv.gz",
    "unique_sequence_map.tsv.gz",
    "rfu_results_per_sequence.tsv.gz",
    "rfu_results_per_row.tsv.gz",
    "run_manifest.json",
)


def _boolean(values: pd.Series) -> pd.Series:
    text = values.astype("string").str.strip().str.lower()
    if (~(text.isin(["true", "false"]) | text.isna())).any():
        raise ValueError("threshold columns contain values other than true/false/missing")
    return text.eq("true").astype("boolean")


def _equal_text(left: pd.Series, right: pd.Series) -> pd.Series:
    return left.astype("string").fillna("<missing>").eq(right.astype("string").fillna("<missing>"))


def validate_completed_rfu_run(
    run_dir: str | Path,
    *,
    evidence_manifest: str | Path | None = None,
    verify_evidence_inputs: bool = False,
    expected_provenance: dict[str, Any] | None = None,
    score_tolerance: float = 1e-12,
    strict: bool = True,
) -> dict[str, Any]:
    """Validate final RFU tables, chunk state, provenance, and optional hashes.

    The validator reads existing tables and manifests only. It never invokes R
    or modifies RFU assignments. CDR3 and chain are sequence-identity features;
    V-call equality is deliberately not required because canonical RFU identity
    is exact CDR3.
    """
    root = Path(run_dir).expanduser().resolve()
    errors: list[str] = []
    warnings: list[str] = []
    missing = [name for name in _REQUIRED_FILES if not (root / name).is_file()]
    if missing:
        errors.append(f"required run files are missing: {missing}")
    if errors:
        report = {"status": "invalid", "errors": errors, "warnings": warnings}
        if strict:
            raise ValueError("Invalid completed RFU run: " + "; ".join(errors))
        return report

    manifest = json.loads((root / "run_manifest.json").read_text())
    receptors = pd.read_csv(root / "receptors.tsv.gz", sep="\t")
    mapping = pd.read_csv(root / "unique_sequence_map.tsv.gz", sep="\t")
    sequences = pd.read_csv(root / "rfu_results_per_sequence.tsv.gz", sep="\t")
    rows = pd.read_csv(root / "rfu_results_per_row.tsv.gz", sep="\t")

    expected_rows = int(manifest.get("original_row_count", -1))
    if not (len(receptors) == len(mapping) == len(rows) == expected_rows):
        errors.append("receptor, mapping, reconstructed, and manifest row counts are inconsistent")
    for label, frame in (("receptors", receptors), ("mapping", mapping), ("rows", rows)):
        if "input_row_id" not in frame or frame["input_row_id"].isna().any():
            errors.append(f"{label} lacks non-missing input_row_id")
        elif frame["input_row_id"].duplicated().any():
            errors.append(f"{label} contains duplicate input_row_id")
    if all("input_row_id" in frame for frame in (receptors, mapping, rows)):
        expected_order = receptors["input_row_id"].astype(str).tolist()
        if mapping["input_row_id"].astype(str).tolist() != expected_order:
            errors.append("unique-sequence mapping does not preserve input order")
        if rows["input_row_id"].astype(str).tolist() != expected_order:
            errors.append("reconstructed output does not preserve input order")

    required_sequence = {
        "unique_sequence_id",
        "cdr3aa",
        "chain",
        "rfu_id",
        "rfu_label",
        "rfu_score",
        "pass_thr",
    }
    if not required_sequence.issubset(sequences):
        errors.append(
            f"per-sequence table is missing columns: {sorted(required_sequence - set(sequences))}"
        )
    elif (
        sequences["unique_sequence_id"].isna().any()
        or sequences["unique_sequence_id"].duplicated().any()
    ):
        errors.append("per-sequence IDs must be unique and non-missing")
    else:
        eligible = mapping.loc[mapping["eligibility_status"].eq("eligible")].copy()
        mapped_ids = set(eligible["unique_sequence_id"].dropna().astype(str))
        sequence_ids = set(sequences["unique_sequence_id"].astype(str))
        if mapped_ids != sequence_ids:
            errors.append("eligible mapping IDs and per-sequence IDs differ")
        intrinsic = (
            eligible.groupby("unique_sequence_id", observed=True)
            .agg(cdr3aa=("cdr3aa", "nunique"))
            .reset_index()
        )
        if intrinsic["cdr3aa"].gt(1).any():
            errors.append("one RFU sequence ID maps to conflicting CDR3 values")
        joined = rows.merge(
            sequences[
                [
                    "unique_sequence_id",
                    "cdr3aa",
                    "chain",
                    "rfu_id",
                    "rfu_label",
                    "rfu_score",
                    "pass_thr",
                ]
            ],
            on="unique_sequence_id",
            how="left",
            suffixes=("_row", "_sequence"),
            sort=False,
            validate="many_to_one",
        )
        eligible_joined = joined.loc[joined["eligibility_status"].eq("eligible")]
        for column in ("cdr3aa", "chain", "rfu_id", "rfu_label"):
            if not _equal_text(
                eligible_joined[f"{column}_row"], eligible_joined[f"{column}_sequence"]
            ).all():
                errors.append(f"per-row and per-sequence {column} values differ")
        row_score = pd.to_numeric(eligible_joined["rfu_score_row"], errors="coerce")
        sequence_score = pd.to_numeric(eligible_joined["rfu_score_sequence"], errors="coerce")
        if not np.allclose(
            row_score.to_numpy(float),
            sequence_score.to_numpy(float),
            rtol=0,
            atol=score_tolerance,
            equal_nan=True,
        ):
            errors.append("per-row and per-sequence RFU scores differ")
        try:
            if not _boolean(eligible_joined["pass_thr_row"]).equals(
                _boolean(eligible_joined["pass_thr_sequence"])
            ):
                errors.append("per-row and per-sequence threshold flags differ")
        except ValueError as error:
            errors.append(str(error))

    if {"rfu_score", "pass_thr", "rfu_id", "assignment_status"}.issubset(rows):
        assigned = rows["rfu_id"].notna()
        score = pd.to_numeric(rows["rfu_score"], errors="coerce")
        try:
            passed = _boolean(rows["pass_thr"]).fillna(False)
            if (passed & (~assigned | score.isna())).any():
                errors.append("threshold pass occurs without a valid nearest assignment and score")
            threshold = float(manifest.get("rfu_threshold", 0.6))
            if (assigned & passed & score.lt(threshold)).any():
                errors.append("threshold-pass rows contain scores below the declared threshold")
            expected_status = pd.Series("ineligible_sequence", index=rows.index, dtype="string")
            eligible = rows["eligibility_status"].eq("eligible")
            expected_status.loc[rows["cdr3aa"].isna()] = "missing_sequence"
            expected_status.loc[eligible & ~assigned] = "upstream_unassigned"
            expected_status.loc[eligible & assigned & passed] = "nearest_threshold_qualified"
            expected_status.loc[eligible & assigned & ~passed] = "nearest_below_threshold"
            if not _equal_text(rows["assignment_status"], expected_status).all():
                errors.append("assignment_status is inconsistent with eligibility/threshold logic")
        except ValueError as error:
            errors.append(str(error))

    backend_manifest_path = Path(str(manifest.get("run_manifest_path", "")))
    if not backend_manifest_path.is_file() and manifest.get("run_id"):
        backend_manifest_path = (
            root / "backend" / "runs" / str(manifest["run_id"]) / "run_manifest.json"
        )
    chunk_count = int(manifest.get("chunk_count", 0))
    completed_chunks = 0
    if chunk_count:
        if not backend_manifest_path.is_file():
            errors.append("backend run manifest is missing")
        else:
            backend_manifest = json.loads(backend_manifest_path.read_text())
            if backend_manifest.get("status") != "complete":
                errors.append("backend run manifest is not complete")
            chunk_dirs = sorted((backend_manifest_path.parent / "chunks").glob("chunk_*"))
            seen_ids: set[str] = set()
            for chunk_dir in chunk_dirs:
                chunk_manifest_path = chunk_dir / "manifest.json"
                output_path = chunk_dir / "output.tsv"
                if not chunk_manifest_path.is_file() or not output_path.is_file():
                    errors.append(f"chunk {chunk_dir.name} is incomplete")
                    continue
                chunk_manifest = json.loads(chunk_manifest_path.read_text())
                chunk_id = str(chunk_manifest.get("chunk_id"))
                if chunk_id in seen_ids:
                    errors.append(f"duplicate chunk ID: {chunk_id}")
                seen_ids.add(chunk_id)
                if chunk_manifest.get("status") != "complete":
                    errors.append(f"chunk {chunk_dir.name} is not complete")
                    continue
                if file_sha256(output_path) != chunk_manifest.get("output_sha256"):
                    errors.append(f"chunk {chunk_dir.name} output checksum mismatch")
                    continue
                completed_chunks += 1
            if len(chunk_dirs) != chunk_count or completed_chunks != chunk_count:
                errors.append("chunk manifest completeness does not match the declared chunk count")

    verified_artifacts = 0
    if evidence_manifest is not None:
        evidence = json.loads(Path(evidence_manifest).expanduser().resolve().read_text())
        groups = ["source_manifests", "outputs"]
        if verify_evidence_inputs:
            groups.append("inputs")
        for group in groups:
            for artifact in evidence.get(group, []):
                path = Path(str(artifact.get("runtime_only_path", "")))
                if not path.is_file():
                    errors.append(f"sealed {group} artifact is missing: {artifact.get('filename')}")
                elif file_sha256(path) != artifact.get("sha256"):
                    errors.append(f"sealed {group} checksum mismatch: {artifact.get('filename')}")
                else:
                    verified_artifacts += 1

    for key, expected in (expected_provenance or {}).items():
        if manifest.get(key) != expected:
            errors.append(f"provenance mismatch for {key!r}")

    try:
        threshold_pass_count = int(_boolean(rows["pass_thr"]).fillna(False).sum())
    except ValueError:
        threshold_pass_count = 0
    report = {
        "status": "valid" if not errors else "invalid",
        "errors": errors,
        "warnings": warnings,
        "run_directory": str(root),
        "row_count": len(rows),
        "eligible_row_count": int(rows["eligibility_status"].eq("eligible").sum()),
        "unique_sequence_count": len(sequences),
        "threshold_pass_count": threshold_pass_count,
        "chunk_count": chunk_count,
        "completed_chunk_count": completed_chunks,
        "verified_evidence_artifact_count": verified_artifacts,
    }
    if strict and errors:
        raise ValueError("Invalid completed RFU run: " + "; ".join(errors))
    return report


__all__ = ["validate_completed_rfu_run"]
