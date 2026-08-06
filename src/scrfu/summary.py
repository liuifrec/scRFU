from __future__ import annotations

import math
import re
from collections.abc import Sequence
from typing import Any

import pandas as pd


def _get_obs_df(adata: Any) -> pd.DataFrame:
    if not hasattr(adata, "obs"):
        raise ValueError("AnnData-like object must have an .obs attribute.")
    return adata.obs


def _require_obs_columns(obs: pd.DataFrame, columns: list[str]) -> None:
    missing = [col for col in columns if col not in obs.columns]
    if missing:
        raise ValueError(f"adata.obs is missing required columns: {missing}")


def _assigned_mask(obs: pd.DataFrame) -> pd.Series:
    labels = obs["rfu_label"]
    return labels.notna()


def _natural_sort_key(value: str) -> list[int | str]:
    parts = re.split(r"(\d+)", str(value))
    return [int(part) if part.isdigit() else part for part in parts]


def _summary_row(df: pd.DataFrame) -> dict[str, object]:
    n_cells = len(df)
    assigned = df.loc[_assigned_mask(df)]
    scores = pd.to_numeric(assigned["rfu_score"], errors="coerce").dropna()
    label_counts = assigned["rfu_label"].value_counts()

    top_rfu = pd.NA
    top_rfu_count = 0
    if not label_counts.empty:
        top_rfu = label_counts.index[0]
        top_rfu_count = int(label_counts.iloc[0])

    return {
        "n_cells": n_cells,
        "n_assigned": len(assigned),
        "assignment_rate": len(assigned) / n_cells if n_cells else 0.0,
        "n_unique_rfu": int(assigned["rfu_label"].nunique(dropna=True)),
        "mean_rfu_score": scores.mean(),
        "median_rfu_score": scores.median(),
        "top_rfu": top_rfu,
        "top_rfu_count": top_rfu_count,
    }


def rfu_summary(adata: Any, groupby: str | None = None) -> pd.DataFrame:
    """
    Summarize RFU assignments stored in ``adata.obs``.

    Uses ``rfu_label`` and ``rfu_score`` columns populated by the RFU calling layer.
    Missing ``rfu_label`` values are treated as unassigned.
    """
    obs = _get_obs_df(adata)
    required = ["rfu_label", "rfu_score"]
    if groupby is None:
        _require_obs_columns(obs, required)
        return pd.DataFrame([_summary_row(obs)])

    _require_obs_columns(obs, [groupby, *required])
    rows = []
    for group_value, group_df in obs.groupby(groupby, dropna=False, sort=True):
        row = {groupby: group_value}
        row.update(_summary_row(group_df))
        rows.append(row)
    return pd.DataFrame(rows)


def aggregate_rfu(adata: Any, groupby: str, normalize: bool = True) -> pd.DataFrame:
    """
    Aggregate assigned RFU labels by group.

    Missing ``rfu_label`` values are excluded from the aggregation.
    """
    obs = _get_obs_df(adata)
    _require_obs_columns(obs, [groupby, "rfu_label"])

    assigned = obs.loc[_assigned_mask(obs), [groupby, "rfu_label"]].copy()
    if assigned.empty:
        return pd.DataFrame(index=pd.Index([], name=groupby))

    table = pd.crosstab(assigned[groupby], assigned["rfu_label"], dropna=True)
    table = table.reindex(sorted(table.columns, key=_natural_sort_key), axis=1)

    if normalize:
        table = table.div(table.sum(axis=1), axis=0)

    table.index.name = groupby
    table.columns.name = None
    return table


def _analysis_frame(data: Any) -> pd.DataFrame:
    if isinstance(data, pd.DataFrame):
        return data.copy()
    if hasattr(data, "obs") and isinstance(data.obs, pd.DataFrame):
        return data.obs.copy()
    raise TypeError("rfu_metrics data must be a pandas DataFrame or an object with DataFrame .obs.")


def _shannon_entropy(counts: pd.Series) -> float:
    total = float(counts.sum())
    if total <= 0:
        return 0.0
    probabilities = counts.astype(float) / total
    return float(-sum(value * math.log(value) for value in probabilities if value > 0))


def rfu_metrics(
    data: Any,
    *,
    groupby: str | Sequence[str],
    weighting: str,
    cell_col: str = "cell_id",
    cdr3_col: str = "cdr3aa",
    rfu_col: str = "rfu_label",
    threshold_col: str = "pass_thr",
    donor_col: str | None = None,
    sample_col: str | None = None,
    chain: str | None = None,
    chain_col: str = "chain",
    assignment_policy: str = "nearest",
) -> pd.DataFrame:
    """Calculate descriptive RFU metrics with explicit phenotype and weighting semantics.

    One output row is produced for each phenotype-group/RFU combination. ``weighting``
    must be ``"cell"`` or ``"unique_sequence"`` and controls ``weighted_abundance``
    and ``threshold_pass_rate``. Counts and the explicitly named cell- and
    sequence-weighted abundance columns are always reported.
    """
    groups = [groupby] if isinstance(groupby, str) else list(groupby)
    if not groups or any(not isinstance(column, str) or not column for column in groups):
        raise ValueError("groupby must name at least one explicit phenotype column.")
    if weighting not in {"cell", "unique_sequence"}:
        raise ValueError("weighting must be explicitly 'cell' or 'unique_sequence'.")
    if assignment_policy not in {"nearest", "threshold_pass"}:
        raise ValueError("assignment_policy must be 'nearest' or 'threshold_pass'.")
    frame = _analysis_frame(data)
    if chain is not None:
        if chain_col not in frame:
            raise ValueError(f"RFU metrics chain selection requires column {chain_col!r}.")
        frame = frame.loc[
            frame[chain_col].astype("string").str.upper().eq(str(chain).upper()).fillna(False)
        ].copy()
    if cell_col not in frame.columns:
        if frame.index.is_unique and frame.index.notna().all():
            frame[cell_col] = frame.index.astype(str)
        else:
            raise ValueError(f"RFU metrics require a {cell_col!r} column or unique index.")
    required = [*groups, cell_col, cdr3_col, rfu_col, threshold_col]
    optional = [column for column in (donor_col, sample_col) if column is not None]
    missing = [column for column in [*required, *optional] if column not in frame.columns]
    if missing:
        raise ValueError(f"RFU metrics input is missing required columns: {missing}")

    assigned = frame[rfu_col].notna() & frame[cdr3_col].notna()
    if assignment_policy == "threshold_pass":
        assigned &= frame[threshold_col].astype("boolean").fillna(False)
    work = frame.loc[assigned, :].copy()
    output_columns = [
        *groups,
        rfu_col,
        "weighting",
        "assignment_policy",
        "rfu_cell_count",
        "rfu_cell_abundance",
        "unique_cdr3_richness",
        "sequence_convergence_ratio",
        "multiplicity",
        "weighted_abundance",
        "clonotype_entropy",
        "dominant_clonotype_fraction",
        "rfu_threshold_pass_rate",
        "cell_abundance",
        "convergence_richness",
        "mean_sequence_multiplicity",
        "normalized_convergence",
        "dominant_sequence_fraction",
        "threshold_pass_rate",
    ]
    for label, column in (("donor", donor_col), ("sample", sample_col)):
        if column is not None:
            output_columns.extend([f"{label}_count", f"{label}_prevalence", f"group_{label}_count"])
    if work.empty:
        return pd.DataFrame(columns=output_columns)

    work[cell_col] = work[cell_col].astype(str)
    work[cdr3_col] = work[cdr3_col].astype(str)
    group_key: str | list[str] = groups[0] if len(groups) == 1 else groups
    metadata_totals: dict[tuple[Any, ...], dict[str, int]] = {}
    for values, group_frame in frame.groupby(group_key, dropna=False, sort=True, observed=True):
        key = values if isinstance(values, tuple) else (values,)
        metadata_totals[key] = {
            f"{label}s": group_frame[column].dropna().nunique()
            for label, column in (("donor", donor_col), ("sample", sample_col))
            if column is not None
        }
    totals: dict[tuple[Any, ...], dict[str, int]] = {}
    for values, group_frame in work.groupby(group_key, dropna=False, sort=True, observed=True):
        key = values if isinstance(values, tuple) else (values,)
        totals[key] = {
            "cells": group_frame[cell_col].nunique(),
            "sequences": group_frame[cdr3_col].nunique(),
            **metadata_totals[key],
        }

    rows: list[dict[str, Any]] = []
    by = [*groups, rfu_col]
    group_rfu_key: str | list[str] = by[0] if len(by) == 1 else by
    for values, subset in work.groupby(group_rfu_key, dropna=False, sort=True, observed=True):
        values_tuple = values if isinstance(values, tuple) else (values,)
        phenotype_values = values_tuple[: len(groups)]
        phenotype_key = tuple(phenotype_values)
        total = totals[phenotype_key]
        observations = subset.drop_duplicates([cell_col, cdr3_col])
        cell_count = subset[cell_col].nunique()
        richness = subset[cdr3_col].nunique()
        observation_count = len(observations)
        clone_counts = observations.groupby(cdr3_col, observed=True)[cell_col].nunique()

        if weighting == "cell":
            weighted_abundance = cell_count / total["cells"] if total["cells"] else 0.0
            pass_units = subset.drop_duplicates(cell_col)
        else:
            weighted_abundance = richness / total["sequences"] if total["sequences"] else 0.0
            pass_units = subset.drop_duplicates(cdr3_col)
        pass_values = pass_units[threshold_col].dropna()
        pass_rate = (
            float(pass_values.astype("boolean").mean()) if not pass_values.empty else float("nan")
        )

        row: dict[str, Any] = {
            **dict(zip(groups, phenotype_values, strict=True)),
            rfu_col: values_tuple[-1],
            "weighting": weighting,
            "assignment_policy": assignment_policy,
            "rfu_cell_count": cell_count,
            "rfu_cell_abundance": cell_count / total["cells"] if total["cells"] else 0.0,
            "unique_cdr3_richness": richness,
            "sequence_convergence_ratio": (
                richness / total["sequences"] if total["sequences"] else 0.0
            ),
            "multiplicity": observation_count / richness if richness else 0.0,
            "weighted_abundance": weighted_abundance,
            "clonotype_entropy": _shannon_entropy(clone_counts),
            "dominant_clonotype_fraction": (
                float(clone_counts.max() / clone_counts.sum()) if not clone_counts.empty else 0.0
            ),
            "rfu_threshold_pass_rate": pass_rate,
            "cell_abundance": cell_count,
            "convergence_richness": richness,
            "mean_sequence_multiplicity": observation_count / richness if richness else 0.0,
            "normalized_convergence": (
                richness / total["sequences"] if total["sequences"] else 0.0
            ),
            "dominant_sequence_fraction": (
                float(clone_counts.max() / clone_counts.sum()) if not clone_counts.empty else 0.0
            ),
            "threshold_pass_rate": pass_rate,
        }
        for label, column in (("donor", donor_col), ("sample", sample_col)):
            if column is not None:
                count = subset[column].dropna().nunique()
                denominator = total[f"{label}s"]
                row[f"{label}_count"] = count
                row[f"{label}_prevalence"] = count / denominator if denominator else float("nan")
                row[f"group_{label}_count"] = denominator
        rows.append(row)
    return pd.DataFrame(rows, columns=output_columns)
