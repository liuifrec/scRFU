from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class StabilityBenchmarkResult:
    metrics: pd.DataFrame
    parameters: dict[str, Any]


def deterministic_subsample(
    data: pd.DataFrame,
    *,
    unit: str = "cell",
    fraction: float | None = None,
    n: int | None = None,
    random_state: int = 0,
    cell_col: str = "cell_id",
    sequence_col: str = "cdr3aa",
    sample_col: str = "sample",
    preserve_order: bool = True,
) -> pd.DataFrame:
    """Subsample whole cells, sequences, or samples without splitting the selected unit."""
    if (fraction is None) == (n is None):
        raise ValueError("Specify exactly one of fraction or n.")
    if fraction is not None and not 0 < fraction <= 1:
        raise ValueError("fraction must lie in (0, 1].")
    if n is not None and (isinstance(n, bool) or not isinstance(n, int) or n < 1):
        raise ValueError("n must be a positive integer.")
    columns = {"cell": cell_col, "sequence": sequence_col, "sample": sample_col}
    if unit not in columns:
        raise ValueError("unit must be 'cell', 'sequence', or 'sample'.")
    column = columns[unit]
    if column not in data:
        raise ValueError(f"Subsampling input is missing {column!r}.")
    units = pd.Index(data[column].dropna().drop_duplicates())
    count = int(np.ceil(len(units) * fraction)) if fraction is not None else min(n or 0, len(units))
    rng = np.random.default_rng(random_state)
    selected = set(rng.choice(units.to_numpy(), size=count, replace=False).tolist())
    result = data.loc[data[column].isin(selected)].copy()
    if not preserve_order:
        result = result.sample(frac=1, random_state=random_state)
    result.attrs["subsampling"] = {
        "unit": unit,
        "selected_unit_count": count,
        "source_unit_count": len(units),
        "fraction": fraction,
        "n": n,
        "random_state": random_state,
        "preserve_order": preserve_order,
    }
    return result


def multinomial_abundance_resample(
    matrix: pd.DataFrame, *, depth: int | pd.Series, random_state: int = 0
) -> pd.DataFrame:
    """Resample abundance vectors; this is not physical read-level downsampling."""
    if matrix.index.has_duplicates or matrix.columns.has_duplicates:
        raise ValueError("Abundance matrix labels must be unique.")
    values = matrix.apply(pd.to_numeric, errors="raise").to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Abundance matrix must contain finite non-negative values.")
    depths = (
        pd.Series(depth, index=matrix.index)
        if np.isscalar(depth)
        else pd.Series(depth).reindex(matrix.index)
    )
    if depths.isna().any() or (depths < 0).any() or not np.allclose(depths, depths.astype(int)):
        raise ValueError("Every resampling depth must be a non-negative integer.")
    rng = np.random.default_rng(random_state)
    output = np.zeros_like(values, dtype=int)
    for index, row in enumerate(values):
        if row.sum() == 0 or depths.iloc[index] == 0:
            continue
        output[index] = rng.multinomial(int(depths.iloc[index]), row / row.sum())
    result = pd.DataFrame(output, index=matrix.index.copy(), columns=matrix.columns.copy())
    result.attrs["resampling"] = {
        "model": "multinomial_abundance",
        "physical_read_downsampling": False,
        "random_state": random_state,
    }
    return result


def _cosine(left: np.ndarray, right: np.ndarray) -> float:
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    return float(np.dot(left, right) / denominator) if denominator else np.nan


def benchmark_representation_stability(
    reference: pd.DataFrame,
    perturbed: pd.DataFrame,
    *,
    metrics: Sequence[str] = (
        "spearman",
        "cosine",
        "jaccard",
        "mean_absolute_error",
        "top_k_overlap",
        "mean_rank_displacement",
    ),
    top_k: int = 10,
) -> StabilityBenchmarkResult:
    """Compare aligned sample representations after one declared perturbation."""
    supported = {
        "spearman",
        "cosine",
        "jaccard",
        "mean_absolute_error",
        "top_k_overlap",
        "mean_rank_displacement",
    }
    invalid = sorted(set(metrics).difference(supported))
    if invalid or top_k < 1:
        raise ValueError(f"Invalid stability metrics {invalid} or top_k.")
    samples = reference.index.intersection(perturbed.index, sort=False)
    features = reference.columns.union(perturbed.columns, sort=False)
    if samples.empty:
        raise ValueError("Representations have no shared samples.")
    left = reference.reindex(index=samples, columns=features, fill_value=0).astype(float)
    right = perturbed.reindex(index=samples, columns=features, fill_value=0).astype(float)
    rows: list[dict[str, Any]] = []
    for sample in samples:
        x, y = left.loc[sample], right.loc[sample]
        rank_x = x.rank(method="average", ascending=False)
        rank_y = y.rank(method="average", ascending=False)
        selected_count = min(top_k, len(features))
        top_x = set(x.sort_values(ascending=False, kind="stable").index[:selected_count])
        top_y = set(y.sort_values(ascending=False, kind="stable").index[:selected_count])
        values = {
            "spearman": float(x.corr(y, method="spearman")),
            "cosine": _cosine(x.to_numpy(), y.to_numpy()),
            "jaccard": len(set(x.index[x.gt(0)]) & set(y.index[y.gt(0)]))
            / len(set(x.index[x.gt(0)]) | set(y.index[y.gt(0)]))
            if (x.gt(0) | y.gt(0)).any()
            else np.nan,
            "mean_absolute_error": float(np.abs(x - y).mean()),
            "top_k_overlap": len(top_x & top_y) / selected_count if selected_count else np.nan,
            "mean_rank_displacement": float(np.abs(rank_x - rank_y).mean()),
        }
        rows.extend(
            {"sample": sample, "metric": metric, "value": values[metric]} for metric in metrics
        )
    return StabilityBenchmarkResult(
        pd.DataFrame(rows),
        {
            "top_k": top_k,
            "shared_sample_count": len(samples),
            "union_feature_count": len(features),
            "metrics": list(metrics),
        },
    )


def threshold_sensitivity(
    data: pd.DataFrame,
    thresholds: Sequence[float],
    *,
    score_col: str = "rfu_score",
    groupby: str | None = None,
) -> pd.DataFrame:
    """Report deterministic threshold-qualified coverage across fixed candidate thresholds."""
    if score_col not in data:
        raise ValueError(f"Threshold sensitivity input is missing {score_col!r}.")
    if groupby is not None and groupby not in data:
        raise ValueError(f"Threshold sensitivity input is missing {groupby!r}.")
    scores = pd.to_numeric(data[score_col], errors="coerce")
    rows: list[dict[str, Any]] = []
    groups = (
        [(None, data)] if groupby is None else data.groupby(groupby, observed=True, dropna=False)
    )
    for threshold in thresholds:
        if not np.isfinite(threshold):
            raise ValueError("Thresholds must be finite.")
        for group, subset in groups:
            eligible_scores = scores.loc[subset.index].dropna()
            row = {
                "threshold": float(threshold),
                "eligible_score_count": len(eligible_scores),
                "threshold_pass_count": int(eligible_scores.ge(threshold).sum()),
                "threshold_pass_fraction": float(eligible_scores.ge(threshold).mean())
                if len(eligible_scores)
                else np.nan,
            }
            if groupby is not None:
                row[groupby] = group
            rows.append(row)
    return pd.DataFrame(rows)


def shuffle_input_order(
    data: pd.DataFrame, *, random_state: int = 0, preserve_index: bool = True
) -> pd.DataFrame:
    """Return a deterministic order perturbation for invariance benchmarks."""
    result = data.sample(frac=1, random_state=random_state)
    if not preserve_index:
        result = result.reset_index(drop=True)
    result.attrs["order_perturbation"] = {
        "random_state": random_state,
        "preserve_index": preserve_index,
    }
    return result


def donor_leave_one_out(data: pd.DataFrame, *, donor_key: str) -> list[tuple[Any, pd.DataFrame]]:
    """Create deterministic donor-held-out folds while preserving donor blocks."""
    if donor_key not in data or data[donor_key].isna().any():
        raise ValueError("Donor leave-one-out requires a non-missing donor column.")
    donors = sorted(data[donor_key].unique(), key=str)
    if len(donors) < 2:
        raise ValueError("Donor leave-one-out requires at least two donors.")
    return [(donor, data.loc[data[donor_key].ne(donor)].copy()) for donor in donors]


__all__ = [
    "StabilityBenchmarkResult",
    "benchmark_representation_stability",
    "deterministic_subsample",
    "donor_leave_one_out",
    "multinomial_abundance_resample",
    "shuffle_input_order",
    "threshold_sensitivity",
]
