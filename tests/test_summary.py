from __future__ import annotations

import pandas as pd
import pytest

from scrfu.tl import aggregate_rfu, rfu_summary


class DummyAdata:
    def __init__(self, obs: pd.DataFrame) -> None:
        self.obs = obs


def test_rfu_summary_global():
    adata = DummyAdata(
        pd.DataFrame(
            {
                "rfu_label": ["RFU10", "RFU2", pd.NA, "RFU10"],
                "rfu_score": [0.1, 0.4, 0.8, 0.7],
                "cluster": ["a", "a", "b", "b"],
            }
        )
    )

    result = rfu_summary(adata)

    assert list(result.columns) == [
        "n_cells",
        "n_assigned",
        "assignment_rate",
        "n_unique_rfu",
        "mean_rfu_score",
        "median_rfu_score",
        "top_rfu",
        "top_rfu_count",
    ]
    assert result.loc[0, "n_cells"] == 4
    assert result.loc[0, "n_assigned"] == 3
    assert result.loc[0, "assignment_rate"] == pytest.approx(0.75)
    assert result.loc[0, "n_unique_rfu"] == 2
    assert result.loc[0, "mean_rfu_score"] == pytest.approx((0.1 + 0.4 + 0.7) / 3)
    assert result.loc[0, "median_rfu_score"] == pytest.approx(0.4)
    assert result.loc[0, "top_rfu"] == "RFU10"
    assert result.loc[0, "top_rfu_count"] == 2


def test_rfu_summary_grouped():
    adata = DummyAdata(
        pd.DataFrame(
            {
                "rfu_label": ["RFU1", pd.NA, "RFU2", "RFU2"],
                "rfu_score": [0.2, 0.5, pd.NA, 0.6],
                "cluster": ["a", "a", "b", "b"],
            }
        )
    )

    result = rfu_summary(adata, groupby="cluster")

    assert list(result["cluster"]) == ["a", "b"]
    assert result.loc[result["cluster"] == "a", "n_assigned"].item() == 1
    assert result.loc[result["cluster"] == "a", "mean_rfu_score"].item() == pytest.approx(0.2)
    assert result.loc[result["cluster"] == "b", "n_assigned"].item() == 2
    assert result.loc[result["cluster"] == "b", "mean_rfu_score"].item() == pytest.approx(0.6)
    assert result.loc[result["cluster"] == "b", "top_rfu"].item() == "RFU2"


def test_aggregate_rfu_normalize_true():
    adata = DummyAdata(
        pd.DataFrame(
            {
                "rfu_label": ["RFU10", "RFU2", pd.NA, "RFU10", "RFU1"],
                "group": ["g1", "g1", "g1", "g2", "g2"],
            }
        )
    )

    result = aggregate_rfu(adata, groupby="group", normalize=True)

    assert list(result.columns) == ["RFU1", "RFU2", "RFU10"]
    assert result.loc["g1"].sum() == pytest.approx(1.0)
    assert result.loc["g2"].sum() == pytest.approx(1.0)
    assert result.loc["g1", "RFU2"] == pytest.approx(0.5)
    assert result.loc["g2", "RFU10"] == pytest.approx(0.5)


def test_aggregate_rfu_normalize_false():
    adata = DummyAdata(
        pd.DataFrame(
            {
                "rfu_label": ["RFU2", "RFU2", "RFU10", pd.NA],
                "group": ["g1", "g1", "g1", "g2"],
            }
        )
    )

    result = aggregate_rfu(adata, groupby="group", normalize=False)

    assert list(result.columns) == ["RFU2", "RFU10"]
    assert result.loc["g1", "RFU2"] == 2
    assert result.loc["g1", "RFU10"] == 1


def test_rfu_summary_requires_rfu_columns():
    adata = DummyAdata(pd.DataFrame({"group": ["a", "b"]}))

    with pytest.raises(ValueError, match="missing required columns"):
        rfu_summary(adata)


def test_aggregate_rfu_requires_groupby_column():
    adata = DummyAdata(pd.DataFrame({"rfu_label": ["RFU1"]}))

    with pytest.raises(ValueError, match="missing required columns"):
        aggregate_rfu(adata, groupby="cluster")


def test_all_cells_unassigned():
    adata = DummyAdata(
        pd.DataFrame(
            {
                "rfu_label": [pd.NA, pd.NA],
                "rfu_score": [0.1, pd.NA],
                "group": ["g1", "g2"],
            }
        )
    )

    summary = rfu_summary(adata)
    aggregate = aggregate_rfu(adata, groupby="group")

    assert summary.loc[0, "n_assigned"] == 0
    assert summary.loc[0, "assignment_rate"] == 0.0
    assert summary.loc[0, "n_unique_rfu"] == 0
    assert pd.isna(summary.loc[0, "top_rfu"])
    assert summary.loc[0, "top_rfu_count"] == 0
    assert aggregate.empty
    assert aggregate.index.name == "group"
