#!/usr/bin/env python3
"""Compare benchmark summaries and warn about large performance regressions."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def compare_baselines(
    baseline: pd.DataFrame,
    current: pd.DataFrame,
    *,
    wall_ratio_warning: float = 1.5,
    rss_ratio_warning: float = 1.5,
) -> pd.DataFrame:
    if wall_ratio_warning <= 0 or rss_ratio_warning <= 0:
        raise ValueError("Warning ratios must be positive.")
    keys = ["operation", "input_rows"]
    required = [*keys, "wall_seconds", "peak_rss_bytes"]
    for label, table in (("baseline", baseline), ("current", current)):
        missing = [column for column in required if column not in table]
        if missing:
            raise ValueError(f"{label} benchmark is missing columns: {missing}")
        if table.duplicated(keys).any():
            raise ValueError(f"{label} benchmark has duplicate operation/input_rows keys.")
    merged = baseline[required].merge(
        current[required], on=keys, how="outer", suffixes=("_baseline", "_current"), indicator=True
    )
    merged["wall_ratio"] = merged["wall_seconds_current"] / merged["wall_seconds_baseline"]
    merged["rss_ratio"] = merged["peak_rss_bytes_current"] / merged["peak_rss_bytes_baseline"]
    merged["status"] = "ok"
    merged.loc[merged["_merge"].eq("left_only"), "status"] = "missing_current"
    merged.loc[merged["_merge"].eq("right_only"), "status"] = "new_operation"
    comparable = merged["_merge"].eq("both")
    merged.loc[
        comparable
        & (merged["wall_ratio"].gt(wall_ratio_warning) | merged["rss_ratio"].gt(rss_ratio_warning)),
        "status",
    ] = "warning"
    return merged.drop(columns="_merge").sort_values(keys, kind="stable").reset_index(drop=True)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--current", type=Path, required=True)
    parser.add_argument("--wall-ratio-warning", type=float, default=1.5)
    parser.add_argument("--rss-ratio-warning", type=float, default=1.5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    result = compare_baselines(
        pd.read_csv(args.baseline, sep="\t"),
        pd.read_csv(args.current, sep="\t"),
        wall_ratio_warning=args.wall_ratio_warning,
        rss_ratio_warning=args.rss_ratio_warning,
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(args.output, sep="\t", index=False)
    print(result.to_string(index=False))
    warning_count = int(result["status"].eq("warning").sum())
    if warning_count:
        print(f"WARNING: {warning_count} benchmark rows exceeded advisory ratios.")


if __name__ == "__main__":
    main()
