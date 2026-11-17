from __future__ import annotations

import itertools
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from .downstream import RFUPseudobulkResult, _pair_metric, rfu_pseudobulk
from .repertoire import analysis_frame


@dataclass(frozen=True)
class LongitudinalDesign:
    design_table: pd.DataFrame
    qc_table: pd.DataFrame
    ordered_donors: tuple[Any, ...]
    ordered_timepoints: tuple[Any, ...]
    parameters: dict[str, Any]
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class RFULongitudinalResult:
    sample_matrix: pd.DataFrame
    donor_time_matrix: pd.DataFrame
    sample_metadata: pd.DataFrame
    qc: pd.DataFrame
    parameters: dict[str, Any]
    missingness_mask: pd.DataFrame


@dataclass(frozen=True)
class LongitudinalDynamicsResult:
    classifications: pd.DataFrame
    trajectories: pd.DataFrame
    missingness_mask: pd.DataFrame
    parameters: dict[str, Any]


@dataclass(frozen=True)
class LongitudinalCompartmentResult:
    pairwise_comparison: pd.DataFrame
    rfu_status: pd.DataFrame
    parameters: dict[str, Any]


@dataclass(frozen=True)
class LongitudinalResamplingResult:
    values: pd.DataFrame
    summary: dict[str, Any]
    parameters: dict[str, Any]


def _time_order(values: pd.Series) -> tuple[list[Any], str]:
    nonmissing = values.dropna()
    if nonmissing.empty:
        raise ValueError("Longitudinal time values are all missing.")
    if isinstance(values.dtype, pd.CategoricalDtype):
        if not values.dtype.ordered:
            raise ValueError("Categorical time values must be explicitly ordered.")
        observed = set(nonmissing)
        return [
            value for value in values.cat.categories if value in observed
        ], "ordered_categorical"
    numeric = pd.to_numeric(nonmissing, errors="coerce")
    if numeric.isna().any():
        raise ValueError(
            "Ambiguous time labels are not coerced. Use numeric values or an ordered categorical."
        )
    lookup = pd.DataFrame({"value": nonmissing, "numeric": numeric}).drop_duplicates()
    conflicts = lookup.groupby("numeric", observed=True)["value"].nunique(dropna=False)
    if (conflicts > 1).any():
        raise ValueError("Distinct time labels map to the same numeric time value.")
    ordered = lookup.sort_values("numeric", kind="stable")["value"].tolist()
    return ordered, "numeric"


def validate_longitudinal_design(
    data: Any,
    *,
    sample_key: str,
    donor_key: str,
    time_key: str,
    compartment_key: str | None = None,
    phenotype_key: str | None = None,
    condition_key: str | None = None,
    tissue_key: str | None = None,
    batch_key: str | None = None,
) -> LongitudinalDesign:
    """Validate a repeated-measures design and return one deterministic row per sample."""
    frame = analysis_frame(data)
    optional = [
        key
        for key in (compartment_key, phenotype_key, condition_key, tissue_key, batch_key)
        if key is not None
    ]
    required = [sample_key, donor_key, time_key, *optional]
    missing = [column for column in required if column not in frame]
    if missing:
        raise ValueError(f"Longitudinal design columns are missing: {missing}")
    if frame[sample_key].isna().any() or frame[donor_key].isna().any():
        raise ValueError("Sample and donor identifiers must be non-missing.")
    if frame[time_key].isna().any():
        raise ValueError(
            "Time values must be non-missing; absent visits belong in the missingness mask."
        )
    if not frame[sample_key].astype(str).map(str.strip).ne("").all():
        raise ValueError("Sample identifiers must be non-empty.")
    ordered_timepoints, time_type = _time_order(frame[time_key])
    mapping_columns = [sample_key, donor_key, time_key, *optional]
    conflicts: list[str] = []
    for column in mapping_columns[1:]:
        counts = frame.groupby(sample_key, observed=True)[column].nunique(dropna=False)
        if (counts > 1).any():
            conflicts.append(column)
    if conflicts:
        raise ValueError(f"Samples map to conflicting longitudinal labels: {conflicts}")
    design = frame[mapping_columns].drop_duplicates(sample_key, keep="first").copy()
    donor_order = sorted(design[donor_key].dropna().unique(), key=str)
    time_rank = {value: index for index, value in enumerate(ordered_timepoints)}
    design["_donor_rank"] = design[donor_key].map({value: i for i, value in enumerate(donor_order)})
    design["_time_rank"] = design[time_key].map(time_rank)
    sort_columns = ["_donor_rank", "_time_rank"]
    if compartment_key:
        sort_columns.append(compartment_key)
    sort_columns.append(sample_key)
    design = design.sort_values(sort_columns, kind="stable").drop(
        columns=["_donor_rank", "_time_rank"]
    )
    qc_rows: list[dict[str, Any]] = []
    warnings: list[str] = []
    for donor in donor_order:
        subset = design.loc[design[donor_key].eq(donor)]
        observed = set(subset[time_key])
        missing_times = [value for value in ordered_timepoints if value not in observed]
        singleton = len(observed) < 2
        if singleton:
            warnings.append(f"Donor {donor!r} has only one observed timepoint.")
        if missing_times:
            warnings.append(f"Donor {donor!r} is missing timepoints: {missing_times}")
        qc_rows.append(
            {
                donor_key: donor,
                "sample_count": int(subset[sample_key].nunique()),
                "observed_timepoint_count": len(observed),
                "missing_timepoint_count": len(missing_times),
                "missing_timepoints": "|".join(map(str, missing_times)),
                "singleton_donor": singleton,
                "repeated_donor": not singleton,
            }
        )
    parameters = {
        "sample_key": sample_key,
        "donor_key": donor_key,
        "time_key": time_key,
        "compartment_key": compartment_key,
        "phenotype_key": phenotype_key,
        "condition_key": condition_key,
        "tissue_key": tissue_key,
        "batch_key": batch_key,
        "time_type": time_type,
    }
    return LongitudinalDesign(
        design.reset_index(drop=True),
        pd.DataFrame(qc_rows),
        tuple(donor_order),
        tuple(ordered_timepoints),
        parameters,
        tuple(warnings),
    )


def _coerce_design(
    design: LongitudinalDesign | pd.DataFrame, parameters: dict[str, Any]
) -> LongitudinalDesign:
    if isinstance(design, LongitudinalDesign):
        return design
    return validate_longitudinal_design(design, **parameters)


def rfu_longitudinal_matrix(
    data: Any,
    *,
    sample_key: str,
    donor_key: str,
    time_key: str,
    design: LongitudinalDesign | pd.DataFrame | None = None,
    compartment_key: str | None = None,
    assignment_policy: str = "nearest",
    weighting: str = "cell",
    normalize: str = "count",
    min_prevalence: int | float = 1,
    pseudocount: float = 0.5,
) -> RFULongitudinalResult:
    """Build sample and donor/time RFU representations without imputing visits."""
    design_parameters = {
        "sample_key": sample_key,
        "donor_key": donor_key,
        "time_key": time_key,
        "compartment_key": compartment_key,
    }
    if design is None:
        if isinstance(data, RFUPseudobulkResult):
            metadata = data.sample_metadata.reset_index()
            design_result = validate_longitudinal_design(metadata, **design_parameters)
        else:
            design_result = validate_longitudinal_design(data, **design_parameters)
    else:
        design_result = _coerce_design(design, design_parameters)
    if isinstance(data, RFUPseudobulkResult):
        matrix = data.matrix.copy()
        pseudobulk_qc = data.qc.copy()
        source_parameters = data.parameters.copy()
    else:
        pseudobulk = rfu_pseudobulk(
            data,
            sample_key=sample_key,
            phenotype_keys=[key for key in (donor_key, time_key, compartment_key) if key],
            assignment_policy=assignment_policy,
            weighting=weighting,
            normalize=normalize,
            min_prevalence=min_prevalence,
            pseudocount=pseudocount,
        )
        matrix = pseudobulk.matrix
        pseudobulk_qc = pseudobulk.qc
        source_parameters = pseudobulk.parameters
    metadata = design_result.design_table.set_index(sample_key)
    unknown = sorted(set(matrix.index).difference(metadata.index), key=str)
    if unknown:
        raise ValueError(f"RFU matrix contains samples absent from longitudinal design: {unknown}")
    absent = sorted(set(metadata.index).difference(matrix.index), key=str)
    if absent:
        raise ValueError(
            "Longitudinal design contains samples absent from the RFU matrix; zero-filling "
            f"would silently impute them: {absent}"
        )
    sample_order = design_result.design_table[sample_key].tolist()
    matrix = matrix.reindex(sample_order)
    metadata = metadata.reindex(sample_order)
    index_columns = [donor_key, time_key, *([compartment_key] if compartment_key else [])]
    donor_time = matrix.copy()
    donor_time.index = pd.MultiIndex.from_frame(
        metadata[index_columns].reset_index(drop=True), names=index_columns
    )
    if donor_time.index.duplicated().any():
        raise ValueError("More than one sample maps to the same donor/time/compartment cell.")
    compartments = (
        sorted(metadata[compartment_key].dropna().unique(), key=str) if compartment_key else [None]
    )
    expected = pd.MultiIndex.from_product(
        [design_result.ordered_donors, design_result.ordered_timepoints, compartments],
        names=[donor_key, time_key, compartment_key or "_compartment"],
    )
    observed = donor_time.index
    if not compartment_key:
        observed = pd.MultiIndex.from_tuples(
            [(*value, None) for value in observed], names=expected.names
        )
    missing_series = pd.Series(~expected.isin(observed), index=expected, name="missing")
    missingness = missing_series.unstack([time_key, compartment_key or "_compartment"])
    parameters = {
        **source_parameters,
        **design_parameters,
        "ordered_timepoints": list(design_result.ordered_timepoints),
        "no_imputation": True,
    }
    return RFULongitudinalResult(
        matrix,
        donor_time,
        metadata,
        pseudobulk_qc,
        parameters,
        missingness,
    )


def longitudinal_similarity(
    data: RFULongitudinalResult | pd.DataFrame,
    *,
    metadata: pd.DataFrame | None = None,
    donor_key: str = "donor",
    time_key: str = "time",
    compartment_key: str | None = None,
    metric: str = "cosine",
) -> pd.DataFrame:
    """Return a tidy sample-pair table with repeated-measure context."""
    if isinstance(data, RFULongitudinalResult):
        matrix = data.sample_matrix
        annotations = data.sample_metadata
        donor_key = data.parameters["donor_key"]
        time_key = data.parameters["time_key"]
        compartment_key = data.parameters.get("compartment_key")
        timepoints = data.parameters["ordered_timepoints"]
    else:
        matrix = data.copy()
        if metadata is None:
            raise ValueError("metadata is required with a bare longitudinal matrix.")
        annotations = metadata.copy()
        timepoints, _ = _time_order(annotations[time_key])
    required = [donor_key, time_key, *([compartment_key] if compartment_key else [])]
    missing = [column for column in required if column not in annotations]
    if missing:
        raise ValueError(f"Longitudinal similarity metadata columns are missing: {missing}")
    annotations = annotations.reindex(matrix.index)
    if annotations[required].isna().any().any():
        raise ValueError("Longitudinal similarity metadata are missing for one or more samples.")
    values = matrix.apply(pd.to_numeric, errors="raise").to_numpy(dtype=float)
    time_rank = {value: index for index, value in enumerate(timepoints)}
    rows: list[dict[str, Any]] = []
    for left in range(len(matrix)):
        for right in range(left + 1, len(matrix)):
            try:
                value = _pair_metric(values[left], values[right], metric)
                status = "valid" if np.isfinite(value) else "undefined"
            except ValueError:
                raise
            sample_a, sample_b = matrix.index[left], matrix.index[right]
            meta_a, meta_b = annotations.iloc[left], annotations.iloc[right]
            time_a, time_b = meta_a[time_key], meta_b[time_key]
            if all(isinstance(value, (int, float, np.number)) for value in (time_a, time_b)):
                interval = abs(float(time_b) - float(time_a))
            else:
                interval = abs(time_rank[time_b] - time_rank[time_a])
            row = {
                "sample_a": sample_a,
                "sample_b": sample_b,
                "donor_a": meta_a[donor_key],
                "donor_b": meta_b[donor_key],
                "time_a": time_a,
                "time_b": time_b,
                "same_donor": meta_a[donor_key] == meta_b[donor_key],
                "time_interval": interval,
                "value": value,
                "metric": metric,
                "direction": "distance"
                if metric in {"bray_curtis_dissimilarity", "jensen_shannon_distance"}
                else "similarity",
                "status": status,
            }
            if compartment_key:
                row.update(
                    {
                        "compartment_a": meta_a[compartment_key],
                        "compartment_b": meta_b[compartment_key],
                        "same_compartment": meta_a[compartment_key] == meta_b[compartment_key],
                    }
                )
            rows.append(row)
    return pd.DataFrame(rows)


def summarize_longitudinal_similarity(
    pairs: pd.DataFrame, *, interval_bins: Sequence[float] | None = None
) -> pd.DataFrame:
    """Summarize valid pair values by donor/compartment relation and optional intervals."""
    required = ["same_donor", "value", "status"]
    missing = [column for column in required if column not in pairs]
    if missing:
        raise ValueError(f"Longitudinal pair table is missing columns: {missing}")
    data = pairs.copy()
    data["donor_relation"] = np.where(data["same_donor"], "within_donor", "between_donor")
    groups = ["donor_relation"]
    if "same_compartment" in data:
        data["compartment_relation"] = np.where(
            data["same_compartment"], "same_compartment", "cross_compartment"
        )
        groups.append("compartment_relation")
    if interval_bins is not None:
        data["time_interval_bin"] = pd.cut(
            data["time_interval"], bins=list(interval_bins), include_lowest=True
        )
        groups.append("time_interval_bin")
    return (
        data.groupby(groups, dropna=False, observed=True)
        .agg(
            pair_count=("value", "size"),
            valid_pair_count=("status", lambda values: int(values.eq("valid").sum())),
            mean_value=("value", "mean"),
            median_value=("value", "median"),
        )
        .reset_index()
    )


def donor_retrieval(
    data: RFULongitudinalResult | pd.DataFrame,
    *,
    metadata: pd.DataFrame | None = None,
    donor_key: str = "donor",
    time_key: str = "time",
    metric: str = "cosine",
    top_k: int = 3,
    exclude_same_timepoint: bool = False,
    compartment_key: str | None = None,
    same_compartment_only: bool = False,
) -> pd.DataFrame:
    """Evaluate donor identity retrieval without fitting parameters on the cohort."""
    if top_k < 1:
        raise ValueError("top_k must be positive.")
    if isinstance(data, RFULongitudinalResult):
        matrix = data.sample_matrix
        annotations = data.sample_metadata
        donor_key = data.parameters["donor_key"]
        time_key = data.parameters["time_key"]
        compartment_key = data.parameters.get("compartment_key") or compartment_key
    else:
        matrix = data.copy()
        if metadata is None:
            raise ValueError("metadata is required with a bare retrieval matrix.")
        annotations = metadata.reindex(matrix.index)
    for column in [donor_key, time_key, *([compartment_key] if compartment_key else [])]:
        if column not in annotations:
            raise ValueError(f"Donor retrieval metadata are missing {column!r}.")
    values = matrix.apply(pd.to_numeric, errors="raise").to_numpy(dtype=float)
    distance = metric in {"bray_curtis_dissimilarity", "jensen_shannon_distance"}
    rows: list[dict[str, Any]] = []
    for query_index, query_sample in enumerate(matrix.index):
        query_meta = annotations.iloc[query_index]
        candidates: list[tuple[int, float]] = []
        for candidate_index in range(len(matrix)):
            if candidate_index == query_index:
                continue
            candidate_meta = annotations.iloc[candidate_index]
            if exclude_same_timepoint and candidate_meta[time_key] == query_meta[time_key]:
                continue
            if (
                same_compartment_only
                and compartment_key
                and (candidate_meta[compartment_key] != query_meta[compartment_key])
            ):
                continue
            value = _pair_metric(values[query_index], values[candidate_index], metric)
            if np.isfinite(value):
                candidates.append((candidate_index, value))
        candidates.sort(
            key=lambda item: (item[1] if distance else -item[1], str(matrix.index[item[0]]))
        )
        correct = [
            rank
            for rank, (index, _) in enumerate(candidates, start=1)
            if annotations.iloc[index][donor_key] == query_meta[donor_key]
        ]
        correct_rank = min(correct) if correct else None
        rows.append(
            {
                "query_sample": query_sample,
                "query_donor": query_meta[donor_key],
                "query_time": query_meta[time_key],
                "candidate_count": len(candidates),
                "correct_donor_rank": correct_rank,
                "top_1_donor_match": correct_rank == 1,
                "top_k_donor_match": correct_rank is not None and correct_rank <= top_k,
                "reciprocal_rank": 1 / correct_rank if correct_rank else 0.0,
                "metric": metric,
                "exclude_same_timepoint": exclude_same_timepoint,
                "same_compartment_only": same_compartment_only,
            }
        )
    result = pd.DataFrame(rows)
    result.attrs["mean_reciprocal_rank"] = (
        float(result["reciprocal_rank"].mean()) if len(result) else np.nan
    )
    return result


def rfu_longitudinal_dynamics(
    data: RFULongitudinalResult,
    *,
    minimum_abundance: float = 0.0,
    minimum_timepoints: int = 2,
    minimum_timepoint_proportion: float = 0.5,
    appearance_threshold: float = 0.0,
    disappearance_threshold: float = 0.0,
    expansion_fold_change: float = 2.0,
    contraction_fold_change: float = 2.0,
    pseudocount: float = 0.5,
) -> LongitudinalDynamicsResult:
    """Classify donor-specific RFU trajectories with explicit descriptive rules."""
    if minimum_timepoints < 1 or not 0 <= minimum_timepoint_proportion <= 1:
        raise ValueError("Invalid minimum longitudinal coverage.")
    if expansion_fold_change <= 1 or contraction_fold_change <= 1 or pseudocount <= 0:
        raise ValueError("Fold-change thresholds must exceed one and pseudocount must be positive.")
    donor_key = data.parameters["donor_key"]
    time_key = data.parameters["time_key"]
    compartment_key = data.parameters.get("compartment_key")
    metadata = data.sample_metadata
    groups = [donor_key, *([compartment_key] if compartment_key else [])]
    trajectory_rows: list[dict[str, Any]] = []
    class_rows: list[dict[str, Any]] = []
    for group_key, sample_meta in metadata.groupby(
        groups[0] if len(groups) == 1 else groups, observed=True, sort=True
    ):
        keys = group_key if isinstance(group_key, tuple) else (group_key,)
        sample_meta = sample_meta.copy()
        rank = {value: index for index, value in enumerate(data.parameters["ordered_timepoints"])}
        sample_meta["_time_rank"] = sample_meta[time_key].map(rank)
        sample_meta = sample_meta.sort_values("_time_rank", kind="stable")
        samples = sample_meta.index.tolist()
        for rfu in data.sample_matrix.columns:
            abundance = data.sample_matrix.loc[samples, rfu].astype(float).to_numpy()
            observed_count = len(abundance)
            coverage = observed_count / len(data.parameters["ordered_timepoints"])
            present = abundance > minimum_abundance
            transitions = (
                int(np.count_nonzero(np.diff(present.astype(int)))) if len(present) > 1 else 0
            )
            first, last = float(abundance[0]), float(abundance[-1])
            fold = (last + pseudocount) / (first + pseudocount)
            sufficiently_observed = (
                observed_count >= minimum_timepoints and coverage >= minimum_timepoint_proportion
            )
            if not sufficiently_observed:
                label = "insufficient_coverage"
            elif first <= appearance_threshold < last:
                label = "appearing"
            elif last <= disappearance_threshold < first:
                label = "disappearing"
            elif fold >= expansion_fold_change:
                label = "expanding"
            elif fold <= 1 / contraction_fold_change:
                label = "contracting"
            elif transitions > 1:
                label = "intermittent"
            elif present.all():
                variability = (abundance.max() + pseudocount) / (abundance.min() + pseudocount)
                label = (
                    "persistent_stable"
                    if variability < max(expansion_fold_change, contraction_fold_change)
                    else "persistent_variable"
                )
            else:
                label = "intermittent"
            identity = dict(zip(groups, keys, strict=True))
            class_rows.append(
                {
                    **identity,
                    "rfu_label": rfu,
                    "classification": label,
                    "observed_timepoint_count": observed_count,
                    "observed_timepoint_proportion": coverage,
                    "first_abundance": first,
                    "last_abundance": last,
                    "fold_change": fold,
                    "presence_transition_count": transitions,
                }
            )
            for sample, (_, meta), value in zip(
                samples, sample_meta.iterrows(), abundance, strict=True
            ):
                trajectory_rows.append(
                    {
                        **identity,
                        "sample": sample,
                        "timepoint": meta[time_key],
                        "rfu_label": rfu,
                        "abundance": value,
                        "observed": True,
                    }
                )
    parameters = {
        "minimum_abundance": minimum_abundance,
        "minimum_timepoints": minimum_timepoints,
        "minimum_timepoint_proportion": minimum_timepoint_proportion,
        "appearance_threshold": appearance_threshold,
        "disappearance_threshold": disappearance_threshold,
        "expansion_fold_change": expansion_fold_change,
        "contraction_fold_change": contraction_fold_change,
        "pseudocount": pseudocount,
        "no_imputation": True,
    }
    return LongitudinalDynamicsResult(
        pd.DataFrame(class_rows),
        pd.DataFrame(trajectory_rows),
        data.missingness_mask.copy(),
        parameters,
    )


def longitudinal_compartment_comparison(
    data: RFULongitudinalResult,
    *,
    metric: str = "weighted_jaccard",
    minimum_abundance: float = 0.0,
) -> LongitudinalCompartmentResult:
    """Compare named compartments within donor/time blocks; names are never hard-coded."""
    compartment_key = data.parameters.get("compartment_key")
    if not compartment_key:
        raise ValueError("Longitudinal compartment comparison requires compartment_key.")
    donor_key = data.parameters["donor_key"]
    time_key = data.parameters["time_key"]
    metadata = data.sample_metadata
    pair_rows: list[dict[str, Any]] = []
    status_rows: list[dict[str, Any]] = []
    for (donor, timepoint), subset in metadata.groupby([donor_key, time_key], observed=True):
        ordered = subset.sort_values(compartment_key, kind="stable")
        for left, right in itertools.combinations(ordered.index, 2):
            left_values = data.sample_matrix.loc[left].to_numpy(dtype=float)
            right_values = data.sample_matrix.loc[right].to_numpy(dtype=float)
            pair_rows.append(
                {
                    donor_key: donor,
                    time_key: timepoint,
                    "sample_a": left,
                    "sample_b": right,
                    "compartment_a": metadata.loc[left, compartment_key],
                    "compartment_b": metadata.loc[right, compartment_key],
                    "metric": metric,
                    "value": _pair_metric(left_values, right_values, metric),
                }
            )
            for rfu, left_value, right_value in zip(
                data.sample_matrix.columns, left_values, right_values, strict=True
            ):
                left_present = left_value > minimum_abundance
                right_present = right_value > minimum_abundance
                status_rows.append(
                    {
                        donor_key: donor,
                        time_key: timepoint,
                        "rfu_label": rfu,
                        "compartment_a": metadata.loc[left, compartment_key],
                        "compartment_b": metadata.loc[right, compartment_key],
                        "abundance_a": left_value,
                        "abundance_b": right_value,
                        "paired_difference": left_value - right_value,
                        "compartment_status": "shared"
                        if left_present and right_present
                        else "specific_a"
                        if left_present
                        else "specific_b"
                        if right_present
                        else "absent",
                    }
                )
    return LongitudinalCompartmentResult(
        pd.DataFrame(pair_rows),
        pd.DataFrame(status_rows),
        {"metric": metric, "minimum_abundance": minimum_abundance, "unit": "donor_time"},
    )


def bootstrap_longitudinal_statistic(
    data: pd.DataFrame,
    statistic: Callable[[pd.DataFrame], float],
    *,
    donor_key: str,
    n_resamples: int = 1000,
    random_state: int = 0,
    confidence_level: float = 0.95,
) -> LongitudinalResamplingResult:
    """Bootstrap donor blocks while preserving every repeated observation in each block."""
    if n_resamples < 1 or not 0 < confidence_level < 1:
        raise ValueError("n_resamples and confidence_level are invalid.")
    if donor_key not in data or data[donor_key].isna().any():
        raise ValueError("Donor bootstrap requires non-missing donor identifiers.")
    donors = sorted(data[donor_key].unique(), key=str)
    if not donors:
        raise ValueError("Donor bootstrap requires at least one donor.")
    rng = np.random.default_rng(random_state)
    values: list[float] = []
    undefined = 0
    for _ in range(n_resamples):
        sampled = rng.choice(donors, size=len(donors), replace=True)
        blocks = []
        for instance, donor in enumerate(sampled):
            block = data.loc[data[donor_key].eq(donor)].copy()
            block["_bootstrap_donor_instance"] = instance
            blocks.append(block)
        value = float(statistic(pd.concat(blocks, ignore_index=True)))
        values.append(value)
        undefined += not np.isfinite(value)
    finite = np.asarray([value for value in values if np.isfinite(value)])
    alpha = (1 - confidence_level) / 2
    summary = {
        "observed": float(statistic(data.copy())),
        "bootstrap_mean": float(finite.mean()) if len(finite) else np.nan,
        "lower_percentile": float(np.quantile(finite, alpha)) if len(finite) else np.nan,
        "upper_percentile": float(np.quantile(finite, 1 - alpha)) if len(finite) else np.nan,
        "undefined_iteration_count": undefined,
    }
    return LongitudinalResamplingResult(
        pd.DataFrame({"iteration": range(n_resamples), "value": values}),
        summary,
        {
            "resampling_unit": "donor",
            "donor_key": donor_key,
            "n_resamples": n_resamples,
            "random_state": random_state,
            "confidence_level": confidence_level,
        },
    )


def permute_longitudinal_labels(
    data: pd.DataFrame,
    statistic: Callable[[pd.DataFrame], float],
    *,
    donor_key: str,
    label_key: str,
    n_permutations: int = 1000,
    random_state: int = 0,
    exact: bool | str = "auto",
) -> LongitudinalResamplingResult:
    """Permute donor-level labels while preserving repeated observations within donors."""
    if n_permutations < 1:
        raise ValueError("n_permutations must be positive.")
    for column in (donor_key, label_key):
        if column not in data:
            raise ValueError(f"Permutation input is missing {column!r}.")
    labels = data.groupby(donor_key, observed=True)[label_key].agg(
        lambda values: values.dropna().unique()
    )
    if any(len(value) != 1 for value in labels):
        raise ValueError("Each donor must map to exactly one non-missing permutation label.")
    donors = sorted(labels.index, key=str)
    observed_labels = [labels.loc[donor][0] for donor in donors]
    use_exact = exact is True or (exact == "auto" and math.factorial(len(donors)) <= n_permutations)
    if exact not in {True, False, "auto"}:
        raise ValueError("exact must be True, False, or 'auto'.")
    if use_exact:
        arrangements = list(dict.fromkeys(itertools.permutations(observed_labels)))
    else:
        rng = np.random.default_rng(random_state)
        arrangements = [tuple(rng.permutation(observed_labels)) for _ in range(n_permutations)]
    values: list[float] = []
    for arrangement in arrangements:
        mapping = dict(zip(donors, arrangement, strict=True))
        permuted = data.copy()
        permuted[label_key] = permuted[donor_key].map(mapping)
        values.append(float(statistic(permuted)))
    finite = np.asarray([value for value in values if np.isfinite(value)])
    observed = float(statistic(data.copy()))
    summary = {
        "observed": observed,
        "null_mean": float(finite.mean()) if len(finite) else np.nan,
        "empirical_upper_tail_probability": (1 + int((finite >= observed).sum()))
        / (1 + len(finite))
        if len(finite)
        else np.nan,
        "undefined_iteration_count": len(values) - len(finite),
    }
    return LongitudinalResamplingResult(
        pd.DataFrame({"iteration": range(len(values)), "value": values}),
        summary,
        {
            "permutation_unit": "donor",
            "donor_key": donor_key,
            "label_key": label_key,
            "requested_permutations": n_permutations,
            "actual_permutations": len(values),
            "random_state": random_state,
            "exact": use_exact,
        },
    )


__all__ = [
    "LongitudinalCompartmentResult",
    "LongitudinalDesign",
    "LongitudinalDynamicsResult",
    "LongitudinalResamplingResult",
    "RFULongitudinalResult",
    "bootstrap_longitudinal_statistic",
    "donor_retrieval",
    "longitudinal_compartment_comparison",
    "longitudinal_similarity",
    "permute_longitudinal_labels",
    "rfu_longitudinal_dynamics",
    "rfu_longitudinal_matrix",
    "summarize_longitudinal_similarity",
    "validate_longitudinal_design",
]
