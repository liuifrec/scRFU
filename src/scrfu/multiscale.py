"""Change retained and hidden by a fixed repertoire grouping.

This is the standard contraction of total variation under aggregation, not a
measure of antigen recognition or functional recovery. Inputs are observed
receptor masses; this module never reconstructs integer counts.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class MultiscaleChangeResult:
    summary: dict[str, Any]
    clone_changes: pd.DataFrame
    group_changes: pd.DataFrame


def _unique_index(values: pd.Series | pd.DataFrame, name: str) -> None:
    if not values.index.is_unique or values.index.hasnans:
        raise ValueError(f"{name} requires unique, non-missing receptor identifiers.")


def multiscale_repertoire_change(
    before: pd.Series,
    after: pd.Series,
    grouping: pd.Series,
    *,
    unmapped: Literal["error", "condition"] = "error",
) -> MultiscaleChangeResult:
    """Compare two observed repertoires using one fixed receptor-to-group map.

    ``before`` and ``after`` contain nonnegative counts or masses indexed by
    the same definition of receptor identity. Absent receptors have zero mass;
    a missing *visit* must not be supplied as an empty observed sample. Values
    are normalized separately to sum one. Empty/zero-mass samples give unknown
    distances, including when both are empty.

    With ``unmapped='condition'``, both distances condition on mapped receptors
    and normalize by their mass. Coverage uses the original totals. Unmapped
    receptors are never pooled into an artificial stable group. The map must
    be fixed across visits; pass a single Series, not visit-specific labels.
    Sampling units and the biological replication unit remain the caller's
    responsibility. Missing mapping entries for zero-mass receptors are ignored.
    """
    if unmapped not in {"error", "condition"}:
        raise ValueError("unmapped must be 'error' or 'condition'.")
    for name, value in (("before", before), ("after", after), ("grouping", grouping)):
        if not isinstance(value, pd.Series):
            raise TypeError(f"{name} must be a pandas Series; missing visits are not samples.")
        _unique_index(value, name)
    index = before.index.union(after.index, sort=False)
    masses = pd.DataFrame(index=index)
    for name, value in (("before", before), ("after", after)):
        numeric = pd.to_numeric(value, errors="raise").astype(float)
        if not np.isfinite(numeric).all() or numeric.lt(0).any():
            raise ValueError("Receptor masses must be finite and nonnegative.")
        masses[name] = numeric.reindex(index, fill_value=0.0)
    masses = masses.loc[masses.sum(axis=1).gt(0)].copy()
    masses["group"] = grouping.reindex(masses.index)
    present = masses["group"].notna()
    if unmapped == "error" and not present.all():
        raise ValueError("Some observed receptors have no group; choose conditioning explicitly.")
    totals = masses[["before", "after"]].sum()
    kept = masses.loc[present].copy()
    assigned = kept[["before", "after"]].sum()
    summary: dict[str, Any] = {
        "status": "valid" if assigned.gt(0).all() else "empty_mapped_sample",
        "estimand": "conditional_on_mapped_receptors"
        if unmapped == "condition"
        else "all_receptors",
        "n_receptors_union": len(masses),
        "n_mapped_receptors_union": len(kept),
        "n_groups_union": int(kept["group"].nunique()),
        "d_clone": float("nan"),
        "d_group": float("nan"),
        "aggregation_cancellation": float("nan"),
    }
    for visit in ("before", "after"):
        summary[f"total_mass_{visit}"] = float(totals[visit])
        summary[f"mapped_mass_{visit}"] = float(assigned[visit])
        summary[f"coverage_{visit}"] = (
            float(assigned[visit] / totals[visit]) if totals[visit] > 0 else float("nan")
        )
    kept["p_before"] = kept["before"] / assigned["before"] if assigned["before"] else np.nan
    kept["p_after"] = kept["after"] / assigned["after"] if assigned["after"] else np.nan
    kept["delta"] = kept["p_after"] - kept["p_before"]
    kept["tv_contribution"] = kept["delta"].abs() / 2
    groups = kept.groupby("group", sort=False, observed=True)[["p_before", "p_after"]].sum(
        min_count=1
    )
    groups["delta"] = groups["p_after"] - groups["p_before"]
    groups["tv_contribution"] = groups["delta"].abs() / 2
    if summary["status"] == "valid":
        d_clone = float(kept["tv_contribution"].sum())
        d_group = float(groups["tv_contribution"].sum())
        if d_group > d_clone + 1e-12:
            raise AssertionError("Fixed grouping must contract total variation.")
        summary.update(
            d_clone=d_clone,
            d_group=d_group,
            aggregation_cancellation=max(0.0, d_clone - d_group),
        )
    return MultiscaleChangeResult(summary, kept, groups)


def permute_fixed_groups(
    grouping: pd.Series,
    *,
    strata: pd.DataFrame | None = None,
    random_state: int = 0,
) -> pd.Series:
    """Shuffle group labels once within fixed receptor-feature strata.

    Group sizes and group-by-stratum counts over the supplied receptor universe
    are preserved exactly. They need not be preserved in each visit's subset.
    Apply the returned map unchanged to every visit. Singleton/uniform strata
    cannot exchange labels: callers should report the exchangeable fraction.
    This is a descriptive grouping control, not a calibrated biological null.
    Sorting identifiers makes results invariant to input row order. Identifiers
    must therefore be strings, with no dependence on outcome-selected inputs.
    """
    _unique_index(grouping, "grouping")
    if grouping.isna().any() or not all(isinstance(x, str) for x in grouping.index):
        raise ValueError("Grouping needs non-missing groups and string receptor identifiers.")
    labels = grouping.sort_index().copy()
    if strata is None:
        features = pd.DataFrame({"all": 0}, index=labels.index)
    else:
        _unique_index(strata, "strata")
        features = strata.reindex(labels.index)
        if features.empty or features.isna().any().any():
            raise ValueError("All receptors need explicit, non-missing matching strata.")
    rng = np.random.default_rng(random_state)
    for indices in features.groupby(
        list(features.columns), sort=True, observed=True
    ).groups.values():
        labels.loc[indices] = rng.permutation(labels.loc[indices].to_numpy())
    return labels.reindex(grouping.index)
