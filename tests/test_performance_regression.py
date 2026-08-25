from __future__ import annotations

import hashlib
import json
from pathlib import Path

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


def test_versioned_performance_baseline_manifest_matches_table() -> None:
    root = Path(__file__).parents[1]
    manifest = json.loads((root / "benchmarks" / "synthetic_python_baseline.json").read_text())
    table_path = root / "benchmarks" / manifest["baseline_file"]
    assert manifest["schema_version"] == 1
    assert hashlib.sha256(table_path.read_bytes()).hexdigest() == manifest["baseline_sha256"]
    table = pd.read_csv(table_path, sep="\t")
    assert sorted(table["input_rows"].unique()) == manifest["sizes"]
    assert set(table["operation"]) == set(manifest["operations"])
    assert manifest["ci_enforcement"] is False
