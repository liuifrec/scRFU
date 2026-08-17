from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd


def analysis_frame(data: Any) -> pd.DataFrame:
    """Return a defensive tabular view of receptor/RFU data."""
    if isinstance(data, pd.DataFrame):
        return data.copy()
    if hasattr(data, "obs") and isinstance(data.obs, pd.DataFrame):
        frame = data.obs.copy()
        if "cell_id" not in frame and frame.index.is_unique and frame.index.notna().all():
            frame["cell_id"] = frame.index.astype(str)
        return frame
    raise TypeError("data must be a pandas DataFrame or an object with DataFrame .obs.")


def shannon(counts: pd.Series | np.ndarray) -> float:
    values = np.asarray(counts, dtype=float)
    total = float(values.sum())
    if total <= 0:
        return 0.0
    probabilities = values[values > 0] / total
    return float(-(probabilities * np.log(probabilities)).sum())


def gini(counts: pd.Series | np.ndarray) -> float:
    values = np.asarray(counts, dtype=float)
    values = values[np.isfinite(values) & (values >= 0)]
    if len(values) == 0 or float(values.sum()) == 0:
        return 0.0
    values.sort()
    n = len(values)
    return float((2 * np.dot(np.arange(1, n + 1), values) / (n * values.sum())) - (n + 1) / n)


def _counts_by(
    subset: pd.DataFrame,
    identity: str,
    *,
    weighting: str,
    cell_col: str,
    cdr3_col: str,
    clonotype_col: str,
) -> pd.Series:
    valid = subset.loc[subset[identity].notna()].copy()
    if valid.empty:
        return pd.Series(dtype=float)
    if weighting == "cell":
        valid = valid.drop_duplicates([cell_col, identity])
        return valid.groupby(identity, observed=True)[cell_col].nunique().astype(float)
    if weighting == "unique_sequence":
        valid = valid.drop_duplicates([cdr3_col, identity])
        return valid.groupby(identity, observed=True)[cdr3_col].nunique().astype(float)
    if clonotype_col not in valid or valid[clonotype_col].isna().all():
        valid = valid.drop_duplicates(cdr3_col)
        return valid.groupby(identity, observed=True)[cdr3_col].nunique().astype(float)
    valid = valid.drop_duplicates([clonotype_col, identity])
    return valid.groupby(identity, observed=True)[clonotype_col].nunique().astype(float)


def repertoire_metrics(
    data: Any,
    *,
    groupby: str | Sequence[str],
    weighting: str,
    cell_col: str = "cell_id",
    cdr3_col: str = "cdr3aa",
    clonotype_col: str = "clonotype_id",
    v_col: str = "v_call",
    productive_col: str = "productive",
    chain: str | None = None,
    chain_col: str = "chain",
) -> pd.DataFrame:
    """Calculate descriptive conventional repertoire metrics by explicit group.

    ``weighting`` is one of ``"cell"``, ``"unique_sequence"``, or
    ``"clonotype"``. A true clonotype column is never synthesized: when it is
    absent, sequence-based fallback fields are named and the fallback is
    recorded explicitly.
    """
    groups = [groupby] if isinstance(groupby, str) else list(groupby)
    if not groups or any(not isinstance(column, str) or not column for column in groups):
        raise ValueError("groupby must name at least one explicit group column.")
    if weighting not in {"cell", "unique_sequence", "clonotype"}:
        raise ValueError("weighting must be 'cell', 'unique_sequence', or 'clonotype'.")
    frame = analysis_frame(data)
    if cell_col not in frame:
        if frame.index.is_unique and frame.index.notna().all():
            frame[cell_col] = frame.index.astype(str)
        else:
            raise ValueError(f"Repertoire metrics require column {cell_col!r} or a unique index.")
    required = [*groups, cell_col, cdr3_col]
    missing = [column for column in required if column not in frame]
    if missing:
        raise ValueError(f"Repertoire metrics input is missing required columns: {missing}")
    if chain is not None:
        if chain_col not in frame:
            raise ValueError(f"Chain selection requires column {chain_col!r}.")
        frame = frame.loc[
            frame[chain_col].astype("string").str.upper().eq(str(chain).upper()).fillna(False)
        ].copy()

    columns = [
        *groups,
        "weighting",
        "diversity_identity",
        "clonotype_fallback",
        "receptor_bearing_row_count",
        "unique_cell_count",
        "unique_cdr3_richness",
        "clonotype_richness",
        "sequence_based_clonotype_richness",
        "shannon_entropy",
        "simpson_diversity",
        "inverse_simpson_diversity",
        "dominant_cdr3_fraction",
        "dominant_clonotype_fraction",
        "mean_cdr3aa_length",
        "median_cdr3aa_length",
        "productive_fraction",
        "trbv_richness",
        "top_trbv_fraction",
        "gini_clonality",
    ]
    if frame.empty:
        return pd.DataFrame(columns=columns)

    key: str | list[str] = groups[0] if len(groups) == 1 else groups
    rows: list[dict[str, Any]] = []
    for values, all_rows in frame.groupby(key, dropna=False, sort=True, observed=True):
        values_tuple = values if isinstance(values, tuple) else (values,)
        subset = all_rows.loc[all_rows[cdr3_col].notna()].copy()
        has_clonotypes = clonotype_col in subset and subset[clonotype_col].notna().any()
        identity = clonotype_col if has_clonotypes else cdr3_col
        clone_counts = _counts_by(
            subset,
            identity,
            weighting=weighting,
            cell_col=cell_col,
            cdr3_col=cdr3_col,
            clonotype_col=clonotype_col,
        )
        cdr3_counts = _counts_by(
            subset,
            cdr3_col,
            weighting=weighting,
            cell_col=cell_col,
            cdr3_col=cdr3_col,
            clonotype_col=clonotype_col,
        )
        total = float(clone_counts.sum())
        sum_squared = float(((clone_counts / total) ** 2).sum()) if total else 0.0
        cdr3_total = float(cdr3_counts.sum())
        if v_col in subset:
            v_counts = _counts_by(
                subset,
                v_col,
                weighting=weighting,
                cell_col=cell_col,
                cdr3_col=cdr3_col,
                clonotype_col=clonotype_col,
            )
        else:
            v_counts = pd.Series(dtype=float)
        lengths = subset[cdr3_col].astype(str).str.len()
        productive_fraction = float("nan")
        if productive_col in all_rows:
            productive = all_rows[productive_col].astype("boolean").dropna()
            if not productive.empty:
                productive_fraction = float(productive.mean())
        rows.append(
            {
                **dict(zip(groups, values_tuple, strict=True)),
                "weighting": weighting,
                "diversity_identity": "clonotype_id" if has_clonotypes else "cdr3aa_sequence",
                "clonotype_fallback": not has_clonotypes,
                "receptor_bearing_row_count": len(subset),
                "unique_cell_count": int(subset[cell_col].nunique()),
                "unique_cdr3_richness": int(subset[cdr3_col].nunique()),
                "clonotype_richness": int(subset[clonotype_col].nunique())
                if has_clonotypes
                else pd.NA,
                "sequence_based_clonotype_richness": int(subset[cdr3_col].nunique()),
                "shannon_entropy": shannon(clone_counts),
                "simpson_diversity": 1.0 - sum_squared if total else 0.0,
                "inverse_simpson_diversity": 1.0 / sum_squared if sum_squared else 0.0,
                "dominant_cdr3_fraction": float(cdr3_counts.max() / cdr3_total)
                if cdr3_total
                else 0.0,
                "dominant_clonotype_fraction": float(clone_counts.max() / total)
                if has_clonotypes and total
                else pd.NA,
                "mean_cdr3aa_length": float(lengths.mean()) if len(lengths) else float("nan"),
                "median_cdr3aa_length": float(lengths.median()) if len(lengths) else float("nan"),
                "productive_fraction": productive_fraction,
                "trbv_richness": int(len(v_counts)),
                "top_trbv_fraction": float(v_counts.max() / v_counts.sum())
                if len(v_counts) and v_counts.sum()
                else float("nan"),
                "gini_clonality": gini(clone_counts),
            }
        )
    return pd.DataFrame(rows, columns=columns)


__all__ = ["analysis_frame", "gini", "repertoire_metrics", "shannon"]
