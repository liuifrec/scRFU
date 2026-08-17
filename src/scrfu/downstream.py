from __future__ import annotations

import math
import re
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from .repertoire import analysis_frame, shannon


def _sort_key(value: Any) -> list[tuple[int, Any]]:
    return [
        (0, int(part)) if part.isdigit() else (1, part.lower())
        for part in re.split(r"(\d+)", str(value))
    ]


def _ordered(values: pd.Series | pd.Index | Sequence[Any]) -> list[Any]:
    unique = pd.Index(values).dropna().drop_duplicates().tolist()
    return sorted(unique, key=_sort_key)


def _require_columns(frame: pd.DataFrame, columns: Sequence[str], context: str) -> None:
    missing = [column for column in columns if column not in frame]
    if missing:
        raise ValueError(f"{context} input is missing required columns: {missing}")


def _assignment_rows(
    frame: pd.DataFrame,
    *,
    assignment_policy: str,
    rfu_col: str,
    threshold_col: str,
) -> pd.DataFrame:
    if assignment_policy not in {"nearest", "threshold_pass"}:
        raise ValueError("assignment_policy must be 'nearest' or 'threshold_pass'.")
    _require_columns(frame, [rfu_col], "RFU assignment")
    mask = frame[rfu_col].notna()
    if assignment_policy == "threshold_pass":
        _require_columns(frame, [threshold_col], "Threshold-pass assignment")
        try:
            mask &= frame[threshold_col].astype("boolean").fillna(False)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{threshold_col!r} must contain boolean-compatible values.") from exc
    return frame.loc[mask].copy()


def _sample_metadata(
    frame: pd.DataFrame, sample_key: str, phenotype_keys: Sequence[str]
) -> pd.DataFrame:
    _require_columns(frame, [sample_key, *phenotype_keys], "RFU pseudobulk")
    if frame[sample_key].isna().any():
        raise ValueError(f"sample_key {sample_key!r} contains missing values.")
    rows: list[dict[str, Any]] = []
    for sample in _ordered(frame[sample_key]):
        subset = frame.loc[frame[sample_key].eq(sample)]
        row: dict[str, Any] = {sample_key: sample}
        for column in phenotype_keys:
            values = subset[column].dropna().drop_duplicates().tolist()
            if len(values) > 1:
                raise ValueError(
                    f"Sample {sample!r} has conflicting {column!r} labels: "
                    f"{sorted(map(str, values))}"
                )
            row[column] = values[0] if values else pd.NA
        rows.append(row)
    return pd.DataFrame(rows).set_index(sample_key)


@dataclass(frozen=True)
class RFUPseudobulkResult:
    matrix: pd.DataFrame
    sample_metadata: pd.DataFrame
    rfu_metadata: pd.DataFrame
    qc: pd.DataFrame
    parameters: dict[str, Any]


@dataclass(frozen=True)
class RFUOverlapResult:
    matrix: pd.DataFrame
    metric: str
    direction: str
    parameters: dict[str, Any]
    sample_metadata: pd.DataFrame | None = None


def rfu_pseudobulk(
    data: Any,
    *,
    sample_key: str,
    phenotype_keys: str | Sequence[str] | None = None,
    assignment_policy: str = "nearest",
    weighting: str = "cell",
    normalize: str = "count",
    min_cells: int = 1,
    min_sequences: int = 1,
    min_prevalence: int | float = 1,
    pseudocount: float = 0.5,
    retain_zero_rfus: bool = False,
    cell_col: str = "cell_id",
    cdr3_col: str = "cdr3aa",
    rfu_col: str = "rfu_label",
    threshold_col: str = "pass_thr",
) -> RFUPseudobulkResult:
    """Build a deterministic biological-sample by RFU matrix."""
    if not isinstance(sample_key, str) or not sample_key:
        raise ValueError("sample_key must explicitly name a biological sample column.")
    phenotypes = (
        []
        if phenotype_keys is None
        else [phenotype_keys]
        if isinstance(phenotype_keys, str)
        else list(phenotype_keys)
    )
    if weighting not in {"cell", "unique_sequence"}:
        raise ValueError("weighting must be 'cell' or 'unique_sequence'.")
    if normalize not in {"count", "proportion", "counts_per_1000", "clr"}:
        raise ValueError("normalize must be 'count', 'proportion', 'counts_per_1000', or 'clr'.")
    if min_cells < 0 or min_sequences < 0:
        raise ValueError("min_cells and min_sequences must be non-negative.")
    if pseudocount <= 0:
        raise ValueError("pseudocount must be positive.")
    frame = analysis_frame(data)
    _require_columns(
        frame, [sample_key, cell_col, cdr3_col, rfu_col, *phenotypes], "RFU pseudobulk"
    )
    metadata = _sample_metadata(frame, sample_key, phenotypes)
    samples = metadata.index.tolist()
    assigned = _assignment_rows(
        frame,
        assignment_policy=assignment_policy,
        rfu_col=rfu_col,
        threshold_col=threshold_col,
    )
    assigned = assigned.loc[assigned[cdr3_col].notna() & assigned[cell_col].notna()].copy()

    qc_rows: list[dict[str, Any]] = []
    eligible_samples: set[Any] = set()
    for sample in samples:
        source_subset = frame.loc[frame[sample_key].eq(sample)]
        subset = assigned.loc[assigned[sample_key].eq(sample)]
        cells = int(subset[cell_col].nunique())
        sequences = int(subset[cdr3_col].nunique())
        passes = cells >= min_cells and sequences >= min_sequences
        if passes:
            eligible_samples.add(sample)
        qc_rows.append(
            {
                sample_key: sample,
                "input_row_count": len(source_subset),
                "assigned_row_count": len(subset),
                "unique_cell_count": cells,
                "unique_sequence_count": sequences,
                "passes_minimum": passes,
            }
        )
    qc = pd.DataFrame(qc_rows).set_index(sample_key)
    assigned = assigned.loc[assigned[sample_key].isin(eligible_samples)]
    unit_col = cell_col if weighting == "cell" else cdr3_col
    units = assigned.drop_duplicates([sample_key, rfu_col, unit_col])

    observed_rfus = _ordered(units[rfu_col]) if len(units) else []
    all_rfus = observed_rfus
    if retain_zero_rfus and isinstance(frame[rfu_col].dtype, pd.CategoricalDtype):
        all_rfus = _ordered([*frame[rfu_col].cat.categories, *observed_rfus])
    counts = pd.DataFrame(0, index=pd.Index(samples, name=sample_key), columns=all_rfus, dtype=int)
    counts.columns.name = rfu_col
    if len(units):
        observed = units.groupby([sample_key, rfu_col], observed=True).size().unstack(fill_value=0)
        observed = observed.reindex(index=samples, columns=all_rfus, fill_value=0)
        counts.loc[:, :] = observed.to_numpy(dtype=int)

    prevalence_count = (counts > 0).sum(axis=0)
    if isinstance(min_prevalence, bool) or not isinstance(min_prevalence, (int, float)):
        raise TypeError("min_prevalence must be an integer sample count or a fraction.")
    if isinstance(min_prevalence, float) and 0 <= min_prevalence <= 1:
        prevalence_required = math.ceil(min_prevalence * len(samples))
    elif float(min_prevalence).is_integer() and min_prevalence >= 0:
        prevalence_required = int(min_prevalence)
    else:
        raise ValueError("min_prevalence must be non-negative; fractional values must be <= 1.")
    prevalence_mask = prevalence_count >= prevalence_required
    if retain_zero_rfus:
        prevalence_mask |= prevalence_count.eq(0)
    kept = prevalence_count.loc[prevalence_mask].index
    counts = counts.loc[:, kept]
    prevalence_count = prevalence_count.loc[kept]

    if normalize == "count":
        matrix = counts
    elif normalize in {"proportion", "counts_per_1000"}:
        denominator = counts.sum(axis=1).replace(0, np.nan)
        scale = 1000.0 if normalize == "counts_per_1000" else 1.0
        matrix = counts.div(denominator, axis=0).mul(scale).fillna(0.0)
    else:
        logged = np.log(counts.astype(float) + pseudocount)
        matrix = logged.sub(logged.mean(axis=1), axis=0) if counts.shape[1] else logged

    rfu_metadata = pd.DataFrame(
        {
            "total_abundance": counts.sum(axis=0),
            "sample_prevalence_count": prevalence_count,
            "sample_prevalence": prevalence_count / len(samples) if samples else np.nan,
        },
        index=pd.Index(counts.columns, name=rfu_col),
    )
    parameters = {
        "sample_key": sample_key,
        "phenotype_keys": phenotypes,
        "assignment_policy": assignment_policy,
        "weighting": weighting,
        "normalize": normalize,
        "min_cells": min_cells,
        "min_sequences": min_sequences,
        "min_prevalence": min_prevalence,
        "pseudocount": pseudocount,
        "retain_zero_rfus": retain_zero_rfus,
        "zero_handling": "additive_pseudocount_before_log" if normalize == "clr" else None,
    }
    return RFUPseudobulkResult(matrix, metadata, rfu_metadata, qc, parameters)


def _pair_metric(left: np.ndarray, right: np.ndarray, metric: str) -> float:
    if np.any(left < 0) or np.any(right < 0):
        raise ValueError("RFU overlap metrics require non-negative abundance values.")
    if metric in {"jaccard", "sorensen_dice", "overlap_coefficient"}:
        a = left > 0
        b = right > 0
        intersection = int(np.logical_and(a, b).sum())
        if metric == "jaccard":
            denominator = int(np.logical_or(a, b).sum())
            return intersection / denominator if denominator else float("nan")
        if metric == "sorensen_dice":
            denominator = int(a.sum() + b.sum())
            return 2 * intersection / denominator if denominator else float("nan")
        denominator = min(int(a.sum()), int(b.sum()))
        return intersection / denominator if denominator else float("nan")
    if metric == "cosine":
        denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
        return float(np.dot(left, right) / denominator) if denominator else float("nan")
    if metric == "bray_curtis_similarity":
        denominator = float((left + right).sum())
        return (
            1.0 - float(np.abs(left - right).sum() / denominator) if denominator else float("nan")
        )
    if metric == "bray_curtis_dissimilarity":
        denominator = float((left + right).sum())
        return float(np.abs(left - right).sum() / denominator) if denominator else float("nan")
    if metric == "weighted_jaccard":
        denominator = float(np.maximum(left, right).sum())
        return float(np.minimum(left, right).sum() / denominator) if denominator else float("nan")
    if metric == "jensen_shannon_distance":
        if left.sum() == 0 or right.sum() == 0:
            return float("nan")
        p = left / left.sum()
        q = right / right.sum()
        midpoint = (p + q) / 2
        divergence = 0.5 * sum(
            np.sum(values[values > 0] * np.log(values[values > 0] / midpoint[values > 0]))
            for values in (p, q)
        )
        return float(math.sqrt(max(divergence, 0.0)))
    raise ValueError(f"Unknown RFU overlap metric: {metric!r}.")


def rfu_overlap(
    data_or_matrix: Any,
    *,
    sample_key: str | None = None,
    metric: str = "jaccard",
    assignment_policy: str = "nearest",
    weighting: str = "cell",
    min_abundance: float = 1,
    **pseudobulk_kwargs: Any,
) -> RFUOverlapResult:
    """Calculate a square sample-by-sample RFU similarity or distance matrix."""
    supported = {
        "jaccard",
        "sorensen_dice",
        "overlap_coefficient",
        "cosine",
        "bray_curtis_similarity",
        "bray_curtis_dissimilarity",
        "weighted_jaccard",
        "jensen_shannon_distance",
    }
    if metric not in supported:
        raise ValueError(f"metric must be one of {sorted(supported)}.")
    if min_abundance < 0:
        raise ValueError("min_abundance must be non-negative.")
    sample_metadata: pd.DataFrame | None
    if isinstance(data_or_matrix, RFUPseudobulkResult):
        matrix = data_or_matrix.matrix.copy()
        sample_metadata = data_or_matrix.sample_metadata.copy()
    elif isinstance(data_or_matrix, pd.DataFrame) and sample_key is None:
        matrix = data_or_matrix.copy()
        sample_metadata = None
    else:
        if sample_key is None:
            raise ValueError("sample_key is required when overlap input is not a matrix.")
        pseudobulk = rfu_pseudobulk(
            data_or_matrix,
            sample_key=sample_key,
            assignment_policy=assignment_policy,
            weighting=weighting,
            normalize="count",
            **pseudobulk_kwargs,
        )
        matrix = pseudobulk.matrix
        sample_metadata = pseudobulk.sample_metadata
    if matrix.index.has_duplicates:
        raise ValueError("RFU overlap matrix sample index must be unique.")
    try:
        values = matrix.apply(pd.to_numeric, errors="raise").to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("RFU overlap matrix must contain numeric values.") from exc
    if not np.isfinite(values).all():
        raise ValueError("RFU overlap matrix must contain only finite values.")
    values[values < min_abundance] = 0.0
    n = len(matrix)
    output = np.empty((n, n), dtype=float)
    for i in range(n):
        for j in range(i, n):
            value = _pair_metric(values[i], values[j], metric)
            output[i, j] = output[j, i] = value
    direction = (
        "distance"
        if metric in {"bray_curtis_dissimilarity", "jensen_shannon_distance"}
        else "similarity"
    )
    result = pd.DataFrame(output, index=matrix.index.copy(), columns=matrix.index.copy())
    return RFUOverlapResult(
        result,
        metric,
        direction,
        {
            "metric": metric,
            "direction": direction,
            "assignment_policy": assignment_policy,
            "weighting": weighting,
            "min_abundance": min_abundance,
        },
        sample_metadata,
    )


def rfu_phenotype_coupling(
    data: Any,
    *,
    phenotype_key: str,
    sample_key: str | None = None,
    assignment_policy: str = "nearest",
    weighting: str = "cell",
    min_abundance: int = 1,
    cell_col: str = "cell_id",
    cdr3_col: str = "cdr3aa",
    rfu_col: str = "rfu_label",
    threshold_col: str = "pass_thr",
) -> pd.DataFrame:
    """Describe RFU abundance, phenotype specificity, and independent-sample recurrence."""
    if weighting not in {"cell", "unique_sequence"}:
        raise ValueError("weighting must be 'cell' or 'unique_sequence'.")
    if min_abundance < 1:
        raise ValueError("min_abundance must be at least 1.")
    frame = analysis_frame(data)
    required = [phenotype_key, cell_col, cdr3_col, rfu_col]
    if sample_key is not None:
        required.append(sample_key)
    _require_columns(frame, required, "RFU phenotype coupling")
    phenotypes = _ordered(frame[phenotype_key])
    assigned = _assignment_rows(
        frame,
        assignment_policy=assignment_policy,
        rfu_col=rfu_col,
        threshold_col=threshold_col,
    )
    assigned = assigned.loc[
        assigned[phenotype_key].notna() & assigned[cell_col].notna() & assigned[cdr3_col].notna()
    ].copy()
    unit_col = cell_col if weighting == "cell" else cdr3_col
    units = assigned.drop_duplicates([phenotype_key, rfu_col, unit_col])
    columns = [
        rfu_col,
        phenotype_key,
        "weighting",
        "assignment_policy",
        "total_abundance",
        "phenotype_specific_abundance",
        "phenotype_specific_proportion",
        "phenotype_entropy",
        "normalized_phenotype_entropy",
        "phenotype_specificity",
        "dominant_phenotype",
        "dominant_phenotype_fraction",
        "represented_phenotype_count",
    ]
    if sample_key is not None:
        columns.extend(
            [
                "sample_recurrence_count",
                "sample_prevalence",
                "phenotype_sample_recurrence_count",
                "phenotype_sample_prevalence",
            ]
        )
    if units.empty:
        return pd.DataFrame(columns=columns)
    counts = units.groupby([rfu_col, phenotype_key], observed=True).size().unstack(fill_value=0)
    counts = counts.reindex(columns=phenotypes, fill_value=0)
    all_sample_count = int(frame[sample_key].dropna().nunique()) if sample_key is not None else 0
    rows: list[dict[str, Any]] = []
    for rfu in _ordered(counts.index):
        abundance = counts.loc[rfu]
        total = int(abundance.sum())
        if total < min_abundance:
            continue
        entropy = shannon(abundance)
        normalized_entropy = entropy / math.log(len(phenotypes)) if len(phenotypes) > 1 else 0.0
        maximum = abundance.max()
        dominant = _ordered(abundance.index[abundance.eq(maximum)])[0]
        rfu_assigned = assigned.loc[assigned[rfu_col].eq(rfu)]
        recurrence = int(rfu_assigned[sample_key].dropna().nunique()) if sample_key else 0
        for phenotype in phenotypes:
            phenotype_abundance = int(abundance.loc[phenotype])
            row: dict[str, Any] = {
                rfu_col: rfu,
                phenotype_key: phenotype,
                "weighting": weighting,
                "assignment_policy": assignment_policy,
                "total_abundance": total,
                "phenotype_specific_abundance": phenotype_abundance,
                "phenotype_specific_proportion": phenotype_abundance / total,
                "phenotype_entropy": entropy,
                "normalized_phenotype_entropy": normalized_entropy,
                "phenotype_specificity": 1.0 - normalized_entropy,
                "dominant_phenotype": dominant,
                "dominant_phenotype_fraction": float(maximum / total),
                "represented_phenotype_count": int((abundance > 0).sum()),
            }
            if sample_key is not None:
                phenotype_samples = int(
                    frame.loc[frame[phenotype_key].eq(phenotype), sample_key].dropna().nunique()
                )
                phenotype_recurrence = int(
                    rfu_assigned.loc[rfu_assigned[phenotype_key].eq(phenotype), sample_key]
                    .dropna()
                    .nunique()
                )
                row.update(
                    {
                        "sample_recurrence_count": recurrence,
                        "sample_prevalence": recurrence / all_sample_count
                        if all_sample_count
                        else float("nan"),
                        "phenotype_sample_recurrence_count": phenotype_recurrence,
                        "phenotype_sample_prevalence": phenotype_recurrence / phenotype_samples
                        if phenotype_samples
                        else float("nan"),
                    }
                )
            rows.append(row)
    return pd.DataFrame(rows, columns=columns)


__all__ = [
    "RFUOverlapResult",
    "RFUPseudobulkResult",
    "rfu_overlap",
    "rfu_phenotype_coupling",
    "rfu_pseudobulk",
]
