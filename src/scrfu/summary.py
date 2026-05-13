from __future__ import annotations

import re
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
