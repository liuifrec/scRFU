from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

from .repertoire import analysis_frame

REFERENCE_COVERAGE_STATUSES = (
    "covered",
    "below_threshold",
    "ineligible_sequence",
    "upstream_unassigned",
    "missing_sequence",
)


def _coverage_status(frame: pd.DataFrame) -> pd.Series:
    if "reference_coverage_status" in frame:
        status = frame["reference_coverage_status"].astype("string")
        invalid = sorted(set(status.dropna()).difference(REFERENCE_COVERAGE_STATUSES))
        if invalid:
            raise ValueError(f"Unknown reference coverage statuses: {invalid}")
        return status
    required = ["cdr3aa", "eligibility_status", "rfu_id", "pass_thr"]
    missing = [column for column in required if column not in frame]
    if missing:
        raise ValueError(
            f"Reference coverage requires reference_coverage_status or legacy columns: {missing}"
        )
    eligible = frame["eligibility_status"].eq("eligible")
    assigned = eligible & frame["rfu_id"].notna()
    passed = frame["pass_thr"].astype("boolean").fillna(False)
    status = pd.Series("ineligible_sequence", index=frame.index, dtype="string")
    status.loc[frame["cdr3aa"].isna()] = "missing_sequence"
    status.loc[eligible & ~assigned] = "upstream_unassigned"
    status.loc[assigned & passed] = "covered"
    status.loc[assigned & ~passed] = "below_threshold"
    return status


def reference_coverage(
    data: Any,
    *,
    groupby: str | Sequence[str] | None = None,
    chain: str | None = None,
    chain_col: str = "chain",
    cell_col: str = "cell_id",
    sequence_col: str = "unique_sequence_id",
    cdr3_col: str = "cdr3aa",
    score_col: str = "rfu_score",
) -> pd.DataFrame:
    """Summarize descriptive frozen-reference coverage without claiming OOD detection."""
    frame = analysis_frame(data)
    groups = [] if groupby is None else [groupby] if isinstance(groupby, str) else list(groupby)
    if len(set(groups)) != len(groups) or any(not column for column in groups):
        raise ValueError("groupby must contain unique, non-empty column names.")
    missing = [column for column in groups if column not in frame]
    if missing:
        raise ValueError(f"Reference coverage grouping columns are missing: {missing}")
    if chain is not None:
        if chain_col not in frame:
            raise ValueError(f"Chain selection requires column {chain_col!r}.")
        frame = frame.loc[
            frame[chain_col].astype("string").str.upper().eq(str(chain).strip().upper())
        ].copy()
    frame["_coverage_status"] = _coverage_status(frame)
    columns = [
        *groups,
        "input_rows",
        "eligible_rows",
        "assigned_rows",
        "threshold_pass_rows",
        "below_threshold_rows",
        "ineligible_rows",
        "upstream_unassigned_rows",
        "missing_sequence_rows",
        "score_q05",
        "score_q25",
        "score_median",
        "score_q75",
        "score_q95",
        "threshold_pass_fraction",
        "unique_sequence_coverage",
        "cell_weighted_coverage",
    ]
    if frame.empty:
        return pd.DataFrame(columns=columns)
    iterator = (
        [((), frame)]
        if not groups
        else frame.groupby(
            groups[0] if len(groups) == 1 else groups,
            dropna=False,
            sort=True,
            observed=True,
        )
    )
    rows: list[dict[str, Any]] = []
    for key, subset in iterator:
        keys = key if isinstance(key, tuple) else (key,)
        status = subset["_coverage_status"]
        eligible = status.isin(["covered", "below_threshold", "upstream_unassigned"])
        assigned = status.isin(["covered", "below_threshold"])
        covered = status.eq("covered")
        scores = pd.to_numeric(
            subset.loc[assigned, score_col] if score_col in subset else pd.Series(dtype=float),
            errors="coerce",
        ).dropna()
        quantiles = scores.quantile([0.05, 0.25, 0.5, 0.75, 0.95]) if len(scores) else None
        if sequence_col in subset:
            eligible_sequences = subset.loc[eligible, sequence_col].dropna().astype(str).nunique()
            covered_sequences = subset.loc[covered, sequence_col].dropna().astype(str).nunique()
        elif cdr3_col in subset:
            eligible_sequences = subset.loc[eligible, cdr3_col].dropna().astype(str).nunique()
            covered_sequences = subset.loc[covered, cdr3_col].dropna().astype(str).nunique()
        else:
            eligible_sequences = covered_sequences = 0
        if cell_col in subset:
            eligible_cells = subset.loc[eligible, cell_col].dropna().astype(str).nunique()
            covered_cells = subset.loc[covered, cell_col].dropna().astype(str).nunique()
        else:
            eligible_cells = covered_cells = 0
        row: dict[str, Any] = dict(zip(groups, keys, strict=True))
        row.update(
            {
                "input_rows": len(subset),
                "eligible_rows": int(eligible.sum()),
                "assigned_rows": int(assigned.sum()),
                "threshold_pass_rows": int(covered.sum()),
                "below_threshold_rows": int(status.eq("below_threshold").sum()),
                "ineligible_rows": int(status.eq("ineligible_sequence").sum()),
                "upstream_unassigned_rows": int(status.eq("upstream_unassigned").sum()),
                "missing_sequence_rows": int(status.eq("missing_sequence").sum()),
                "score_q05": float(quantiles.loc[0.05]) if quantiles is not None else np.nan,
                "score_q25": float(quantiles.loc[0.25]) if quantiles is not None else np.nan,
                "score_median": float(quantiles.loc[0.5]) if quantiles is not None else np.nan,
                "score_q75": float(quantiles.loc[0.75]) if quantiles is not None else np.nan,
                "score_q95": float(quantiles.loc[0.95]) if quantiles is not None else np.nan,
                "threshold_pass_fraction": float(covered.sum() / eligible.sum())
                if eligible.any()
                else np.nan,
                "unique_sequence_coverage": covered_sequences / eligible_sequences
                if eligible_sequences
                else np.nan,
                "cell_weighted_coverage": covered_cells / eligible_cells
                if eligible_cells
                else np.nan,
            }
        )
        rows.append(row)
    return pd.DataFrame(rows, columns=columns)


__all__ = ["REFERENCE_COVERAGE_STATUSES", "reference_coverage"]
