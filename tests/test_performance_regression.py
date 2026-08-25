from __future__ import annotations

import pandas as pd
import pytest

from examples.compare_benchmark_baseline import compare_baselines


def test_performance_comparison_is_advisory_and_keyed() -> None:
    baseline = pd.DataFrame(
        {
            "operation": ["validate", "overlap"],
            "input_rows": [100, 100],
            "wall_seconds": [1.0, 2.0],
            "peak_rss_bytes": [1000, 2000],
        }
    )
    current = baseline.copy()
    current.loc[0, "wall_seconds"] = 2.0
    result = compare_baselines(baseline, current)
    assert result.set_index("operation").loc["validate", "status"] == "warning"
    assert result.set_index("operation").loc["overlap", "status"] == "ok"


def test_performance_comparison_rejects_ambiguous_or_incomplete_keys() -> None:
    table = pd.DataFrame(
        {
            "operation": ["x", "x"],
            "input_rows": [10, 10],
            "wall_seconds": [1, 1],
            "peak_rss_bytes": [1, 1],
        }
    )
    with pytest.raises(ValueError, match="duplicate"):
        compare_baselines(table, table)
