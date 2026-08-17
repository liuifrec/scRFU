from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from typing import Any

import numpy as np
import pandas as pd

from ._version import __version__
from .diagnostics import reference_coverage
from .downstream import rfu_pseudobulk
from .repertoire import analysis_frame

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_HARMONIZED_FIELDS = {
    "cohort",
    "sample",
    "donor",
    "age",
    "sex",
    "tissue",
    "cell_subset",
    "condition",
    "platform",
    "sequencing_depth",
    "chain",
    "timepoint",
}


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _canonical_reference_payload(values: Mapping[str, Any]) -> dict[str, Any]:
    return {
        key: values[key]
        for key in (
            "rfu_r_sha256",
            "embedding_sha256",
            "centroid_atlas_sha256",
            "threshold",
            "eligibility_rule",
            "assignment_mode",
            "receptor_chain",
            "receptor_model",
            "reference_label",
        )
    }


def _reference_identifier(values: Mapping[str, Any]) -> str:
    payload = json.dumps(
        _canonical_reference_payload(values), sort_keys=True, separators=(",", ":")
    ).encode()
    return "scrfu-ref-" + hashlib.sha256(payload).hexdigest()


@dataclass(frozen=True)
class FrozenRFUReference:
    rfu_r_sha256: str
    embedding_sha256: str
    centroid_atlas_sha256: str
    threshold: float
    eligibility_rule: str
    assignment_mode: str
    receptor_chain: str
    receptor_model: str
    reference_label: str
    immutable_reference_id: str
    creation_timestamp: str
    scrfu_version: str
    schema_version: str = "1"

    @classmethod
    def create(
        cls,
        *,
        rfu_r_sha256: str,
        embedding_sha256: str,
        centroid_atlas_sha256: str,
        threshold: float,
        eligibility_rule: str,
        assignment_mode: str,
        receptor_chain: str,
        receptor_model: str,
        reference_label: str,
        creation_timestamp: str | None = None,
    ) -> FrozenRFUReference:
        values = {
            "rfu_r_sha256": rfu_r_sha256,
            "embedding_sha256": embedding_sha256,
            "centroid_atlas_sha256": centroid_atlas_sha256,
            "threshold": threshold,
            "eligibility_rule": eligibility_rule,
            "assignment_mode": assignment_mode,
            "receptor_chain": receptor_chain,
            "receptor_model": receptor_model,
            "reference_label": reference_label,
        }
        reference = cls(
            **values,
            immutable_reference_id=_reference_identifier(values),
            creation_timestamp=creation_timestamp or _utc_now(),
            scrfu_version=__version__,
        )
        validate_frozen_reference(reference)
        return reference

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class TransferCohortResult:
    reference: FrozenRFUReference
    coverage: pd.DataFrame
    score_distribution: pd.DataFrame
    rfu_summary: pd.DataFrame
    sample_matrix: pd.DataFrame
    sample_metadata: pd.DataFrame
    source_table: pd.DataFrame
    provenance: dict[str, Any]


@dataclass(frozen=True)
class CohortHarmonizationResult:
    metadata: pd.DataFrame
    qc: pd.DataFrame
    unmapped_values: pd.DataFrame
    rules: dict[str, Any]


@dataclass(frozen=True)
class HeldOutValidationManifest:
    development_cohorts: tuple[str, ...]
    validation_cohorts: tuple[str, ...]
    held_out_cohort: str
    frozen_reference_id: str
    frozen_parameters: dict[str, Any]
    data_hashes: dict[str, str]
    evaluation_metrics: tuple[str, ...]
    code_version: str
    created_at: str
    schema_version: str = "1"

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["development_cohorts"] = list(self.development_cohorts)
        result["validation_cohorts"] = list(self.validation_cohorts)
        result["evaluation_metrics"] = list(self.evaluation_metrics)
        return result


def validate_frozen_reference(
    reference: FrozenRFUReference | Mapping[str, Any],
) -> FrozenRFUReference:
    """Validate artifact hashes and the identifier of an immutable RFU reference."""
    if not isinstance(reference, FrozenRFUReference):
        try:
            reference = FrozenRFUReference(**dict(reference))
        except (TypeError, ValueError) as exc:
            raise ValueError("Invalid frozen RFU reference fields.") from exc
    for field in ("rfu_r_sha256", "embedding_sha256", "centroid_atlas_sha256"):
        value = getattr(reference, field)
        if not isinstance(value, str) or not _SHA256.fullmatch(value.lower()):
            raise ValueError(f"{field} must be a 64-character hexadecimal SHA256.")
    if not np.isfinite(reference.threshold):
        raise ValueError("Frozen-reference threshold must be finite.")
    if not all(
        isinstance(getattr(reference, field), str) and getattr(reference, field).strip()
        for field in (
            "eligibility_rule",
            "assignment_mode",
            "receptor_chain",
            "receptor_model",
            "reference_label",
            "creation_timestamp",
            "scrfu_version",
        )
    ):
        raise ValueError("Frozen-reference text fields must be non-empty.")
    expected = _reference_identifier(asdict(reference))
    if reference.immutable_reference_id != expected:
        raise ValueError("Frozen-reference identifier does not match its immutable parameters.")
    return reference


def transfer_cohort(
    data: Any,
    reference: FrozenRFUReference | Mapping[str, Any],
    *,
    cohort_label: str,
    sample_key: str,
    assignment_policy: str = "nearest",
    weighting: str = "cell",
    groupby: str | Sequence[str] | None = None,
    observed_reference_id: str | None = None,
    cell_col: str = "cell_id",
    cdr3_col: str = "cdr3aa",
    rfu_col: str = "rfu_label",
    chain_col: str = "chain",
) -> TransferCohortResult:
    """Summarize already assigned target-cohort rows against one frozen reference.

    This function never fits centroids, changes a threshold, or executes the RFU backend.
    """
    frozen = validate_frozen_reference(reference)
    if observed_reference_id is not None and observed_reference_id != frozen.immutable_reference_id:
        raise ValueError("Target results were not generated with the supplied frozen reference.")
    if not isinstance(cohort_label, str) or not cohort_label.strip():
        raise ValueError("cohort_label must be explicit and non-empty.")
    frame = analysis_frame(data).copy()
    required = [sample_key, cell_col, cdr3_col, rfu_col]
    missing = [column for column in required if column not in frame]
    if missing:
        raise ValueError(f"Transfer cohort source table is missing columns: {missing}")
    chain_verified = False
    if chain_col in frame:
        observed_chains = sorted(
            frame[chain_col].dropna().astype(str).str.strip().str.upper().unique()
        )
        incompatible = [
            chain for chain in observed_chains if chain != frozen.receptor_chain.strip().upper()
        ]
        if incompatible:
            raise ValueError(
                "Target receptor chains are incompatible with the frozen reference: "
                f"{incompatible} versus {frozen.receptor_chain!r}."
            )
        chain_verified = bool(observed_chains)
    coverage_groups = [sample_key]
    if groupby is not None:
        extras = [groupby] if isinstance(groupby, str) else list(groupby)
        coverage_groups.extend(column for column in extras if column not in coverage_groups)
    coverage = reference_coverage(frame, groupby=coverage_groups)
    pseudobulk = rfu_pseudobulk(
        frame,
        sample_key=sample_key,
        assignment_policy=assignment_policy,
        weighting=weighting,
        normalize="count",
        cell_col=cell_col,
        cdr3_col=cdr3_col,
        rfu_col=rfu_col,
    )
    score = pd.to_numeric(
        frame["rfu_score"] if "rfu_score" in frame else pd.Series(dtype=float),
        errors="coerce",
    ).dropna()
    score_distribution = pd.DataFrame(
        {
            "quantile": [0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0],
            "rfu_score": score.quantile([0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0]).to_numpy()
            if len(score)
            else [np.nan] * 7,
        }
    )
    assigned = frame.loc[frame[rfu_col].notna()]
    rfu_summary = (
        assigned.groupby(rfu_col, observed=True)
        .agg(
            row_count=(cell_col, "size"),
            cell_count=(cell_col, "nunique"),
            unique_sequence_count=(cdr3_col, "nunique"),
            sample_prevalence=(sample_key, "nunique"),
        )
        .reset_index()
        .sort_values(rfu_col, kind="stable")
        .reset_index(drop=True)
    )
    unique_sequences = int(frame[cdr3_col].dropna().astype(str).nunique())
    provenance = {
        "schema_version": "1",
        "cohort_label": cohort_label,
        "frozen_reference_id": frozen.immutable_reference_id,
        "scRFU_version": __version__,
        "assignment_policy": assignment_policy,
        "weighting": weighting,
        "source_row_count": len(frame),
        "unique_sequence_count": unique_sequences,
        "deduplication_ratio": len(frame) / unique_sequences if unique_sequences else np.nan,
        "created_at": _utc_now(),
        "target_reference_verified": observed_reference_id is not None,
        "target_chain_verified": chain_verified,
        "no_reference_refitting": True,
    }
    return TransferCohortResult(
        frozen,
        coverage,
        score_distribution,
        rfu_summary,
        pseudobulk.matrix,
        pseudobulk.sample_metadata,
        frame,
        provenance,
    )


def harmonize_cohort_metadata(
    metadata: pd.DataFrame,
    *,
    cohort_label: str,
    field_mapping: Mapping[str, str],
    value_mappings: Mapping[str, Mapping[Any, Any]] | None = None,
    strict: bool = True,
) -> CohortHarmonizationResult:
    """Apply only explicit source-column and source-value mappings."""
    if not isinstance(metadata, pd.DataFrame):
        raise TypeError("metadata must be a pandas DataFrame.")
    if not cohort_label:
        raise ValueError("cohort_label must be explicit.")
    unknown = sorted(set(field_mapping).difference(_HARMONIZED_FIELDS))
    if unknown:
        raise ValueError(f"Unknown harmonized fields: {unknown}")
    missing = sorted(set(field_mapping.values()).difference(metadata.columns))
    if missing:
        raise ValueError(f"Mapped source columns are missing: {missing}")
    duplicate_sources = pd.Series(list(field_mapping.values())).duplicated(keep=False)
    if duplicate_sources.any():
        duplicates = sorted(pd.Series(list(field_mapping.values()))[duplicate_sources].unique())
        raise ValueError(f"Source columns map to multiple canonical fields: {duplicates}")
    mappings = dict(value_mappings or {})
    unexpected_value_fields = sorted(set(mappings).difference(field_mapping))
    if unexpected_value_fields:
        raise ValueError(f"Value mappings lack field mappings: {unexpected_value_fields}")
    result = pd.DataFrame(index=metadata.index.copy())
    result["cohort"] = cohort_label
    qc_rows: list[dict[str, Any]] = []
    unmapped_rows: list[dict[str, Any]] = []
    for canonical, source in field_mapping.items():
        values = metadata[source]
        mapping = mappings.get(canonical)
        if mapping is None:
            result[canonical] = values
            unmapped_count = 0
        else:
            normalized = values.map(mapping)
            unmapped_mask = values.notna() & normalized.isna()
            for source_value, count in values.loc[unmapped_mask].value_counts(dropna=False).items():
                unmapped_rows.append(
                    {
                        "cohort": cohort_label,
                        "field": canonical,
                        "source_field": source,
                        "source_value": source_value,
                        "row_count": int(count),
                    }
                )
            if strict and unmapped_mask.any():
                examples = values.loc[unmapped_mask].drop_duplicates().astype(str).tolist()[:5]
                raise ValueError(f"Unmapped values for {canonical!r}: {examples}")
            result[canonical] = normalized
            unmapped_count = int(unmapped_mask.sum())
        qc_rows.append(
            {
                "cohort": cohort_label,
                "field": canonical,
                "source_field": source,
                "row_count": len(values),
                "missing_count": int(values.isna().sum()),
                "unmapped_count": unmapped_count,
                "source_unique_count": int(values.nunique(dropna=True)),
                "harmonized_unique_count": int(result[canonical].nunique(dropna=True)),
            }
        )
    source_only = [column for column in metadata if column not in field_mapping.values()]
    rules = {
        "cohort_label": cohort_label,
        "field_mapping": dict(field_mapping),
        "value_mappings": {key: dict(value) for key, value in mappings.items()},
        "source_only_fields": source_only,
        "strict": strict,
        "no_inferred_mappings": True,
    }
    return CohortHarmonizationResult(
        result,
        pd.DataFrame(qc_rows),
        pd.DataFrame(
            unmapped_rows,
            columns=["cohort", "field", "source_field", "source_value", "row_count"],
        ),
        rules,
    )


def create_heldout_validation_manifest(
    *,
    development_cohorts: Sequence[str],
    validation_cohorts: Sequence[str],
    held_out_cohort: str,
    reference: FrozenRFUReference | Mapping[str, Any],
    frozen_parameters: Mapping[str, Any],
    data_hashes: Mapping[str, str],
    evaluation_metrics: Sequence[str],
) -> HeldOutValidationManifest:
    """Freeze cohort roles, inputs, parameters, and metrics before held-out evaluation."""
    frozen = validate_frozen_reference(reference)
    development = tuple(map(str, development_cohorts))
    validation = tuple(map(str, validation_cohorts))
    held_out = str(held_out_cohort)
    if not development or not held_out or not evaluation_metrics:
        raise ValueError("Development cohorts, held-out cohort, and metrics are required.")
    overlap = set(development).intersection(validation) | {held_out}.intersection(
        set(development) | set(validation)
    )
    if overlap:
        raise ValueError(f"Cohort roles overlap: {sorted(overlap)}")
    cohorts = set(development) | set(validation) | {held_out}
    missing_hashes = sorted(cohorts.difference(data_hashes))
    extra_hashes = sorted(set(data_hashes).difference(cohorts))
    if missing_hashes or extra_hashes:
        raise ValueError(
            f"Data hashes must exactly cover registered cohorts; missing={missing_hashes}, "
            f"extra={extra_hashes}"
        )
    invalid_hashes = [
        label for label, value in data_hashes.items() if not _SHA256.fullmatch(str(value).lower())
    ]
    if invalid_hashes:
        raise ValueError(f"Invalid cohort data SHA256 values: {invalid_hashes}")
    return HeldOutValidationManifest(
        development,
        validation,
        held_out,
        frozen.immutable_reference_id,
        dict(frozen_parameters),
        dict(data_hashes),
        tuple(map(str, evaluation_metrics)),
        __version__,
        _utc_now(),
    )


__all__ = [
    "CohortHarmonizationResult",
    "FrozenRFUReference",
    "HeldOutValidationManifest",
    "TransferCohortResult",
    "create_heldout_validation_manifest",
    "harmonize_cohort_metadata",
    "transfer_cohort",
    "validate_frozen_reference",
]
