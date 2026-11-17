from __future__ import annotations

import math
import re
from collections.abc import Sequence
from typing import Any

import numpy as np
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
    by = [*groups, rfu_col]
    grouped = work.groupby(by, dropna=False, sort=True, observed=True)
    metrics = pd.concat(
        [
            grouped[cell_col].nunique().rename("rfu_cell_count"),
            grouped[cdr3_col].nunique().rename("unique_cdr3_richness"),
        ],
        axis=1,
    )

    observations = work.drop_duplicates([*by, cell_col, cdr3_col])
    observation_counts = observations.groupby(by, dropna=False, sort=True, observed=True).size()
    clone_counts = observations.groupby([*by, cdr3_col], dropna=False, sort=True, observed=True)[
        cell_col
    ].nunique()
    group_levels = list(range(len(by)))
    clone_totals = clone_counts.groupby(level=group_levels, dropna=False).transform("sum")
    probabilities = clone_counts.astype(float) / clone_totals
    entropy = (
        (-(probabilities * np.log(probabilities))).groupby(level=group_levels, dropna=False).sum()
    )
    dominant = (
        clone_counts.groupby(level=group_levels, dropna=False).max()
        / clone_counts.groupby(level=group_levels, dropna=False).sum()
    )
    metrics["observation_count"] = observation_counts
    metrics["clonotype_entropy"] = entropy
    metrics["dominant_clonotype_fraction"] = dominant

    pass_unit = cell_col if weighting == "cell" else cdr3_col
    pass_units = work.drop_duplicates([*by, pass_unit]).copy()
    pass_units["_threshold_boolean"] = pass_units[threshold_col].astype("boolean")
    pass_rate = pass_units.groupby(by, dropna=False, sort=True, observed=True)[
        "_threshold_boolean"
    ].mean()
    metrics["rfu_threshold_pass_rate"] = pass_rate
    metrics = metrics.reset_index()

    group_totals = (
        work.groupby(groups, dropna=False, sort=True, observed=True)
        .agg(
            group_cell_count=(cell_col, "nunique"),
            group_sequence_count=(cdr3_col, "nunique"),
        )
        .reset_index()
    )
    metrics = metrics.merge(group_totals, on=groups, how="left", sort=False, validate="many_to_one")

    for label, column in (("donor", donor_col), ("sample", sample_col)):
        if column is None:
            continue
        group_column = f"group_{label}_count"
        count_column = f"{label}_count"
        metadata_totals = (
            frame.groupby(groups, dropna=False, sort=True, observed=True)[column]
            .nunique()
            .rename(group_column)
            .reset_index()
        )
        rfu_counts = (
            work.groupby(by, dropna=False, sort=True, observed=True)[column]
            .nunique()
            .rename(count_column)
            .reset_index()
        )
        metrics = metrics.merge(
            metadata_totals, on=groups, how="left", sort=False, validate="many_to_one"
        ).merge(rfu_counts, on=by, how="left", sort=False, validate="one_to_one")
        metrics[f"{label}_prevalence"] = metrics[count_column] / metrics[group_column].replace(
            0, np.nan
        )

    cell_fraction = metrics["rfu_cell_count"] / metrics["group_cell_count"].replace(0, np.nan)
    sequence_fraction = metrics["unique_cdr3_richness"] / metrics["group_sequence_count"].replace(
        0, np.nan
    )
    multiplicity = metrics["observation_count"] / metrics["unique_cdr3_richness"].replace(0, np.nan)
    metrics["weighting"] = weighting
    metrics["assignment_policy"] = assignment_policy
    metrics["rfu_cell_abundance"] = cell_fraction.fillna(0.0)
    metrics["sequence_convergence_ratio"] = sequence_fraction.fillna(0.0)
    metrics["multiplicity"] = multiplicity.fillna(0.0)
    metrics["weighted_abundance"] = (
        metrics["rfu_cell_abundance"]
        if weighting == "cell"
        else metrics["sequence_convergence_ratio"]
    )
    metrics["cell_abundance"] = metrics["rfu_cell_count"]
    metrics["convergence_richness"] = metrics["unique_cdr3_richness"]
    metrics["mean_sequence_multiplicity"] = metrics["multiplicity"]
    metrics["normalized_convergence"] = metrics["sequence_convergence_ratio"]
    metrics["dominant_sequence_fraction"] = metrics["dominant_clonotype_fraction"]
    metrics["threshold_pass_rate"] = metrics["rfu_threshold_pass_rate"]
    return metrics.loc[:, output_columns]
