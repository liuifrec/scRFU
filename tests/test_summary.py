from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from scrfu.tl import aggregate_rfu, rfu_metrics, rfu_summary


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


def _metrics_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4", "c5", "c6"],
            "phenotype": ["P1", "P1", "P1", "P1", "P1", "P2"],
            "cdr3aa": ["A", "A", "B", "C", "D", "E"],
            "rfu_label": ["RFU1", "RFU1", "RFU1", "RFU2", "RFU2", "RFU1"],
            "pass_thr": [True, True, False, True, False, True],
            "donor": ["d1", "d1", "d2", "d2", "d3", "d4"],
            "sample": ["s1", "s1", "s2", "s2", "s3", "s4"],
        }
    )


def test_rfu_metrics_cell_weighting_and_prevalence() -> None:
    result = rfu_metrics(
        _metrics_frame(),
        groupby="phenotype",
        weighting="cell",
        donor_col="donor",
        sample_col="sample",
    )
    row = result[(result["phenotype"] == "P1") & (result["rfu_label"] == "RFU1")].iloc[0]

    assert row["rfu_cell_count"] == 3
    assert row["rfu_cell_abundance"] == pytest.approx(3 / 5)
    assert row["unique_cdr3_richness"] == 2
    assert row["sequence_convergence_ratio"] == pytest.approx(2 / 4)
    assert row["multiplicity"] == pytest.approx(3 / 2)
    assert row["weighted_abundance"] == pytest.approx(3 / 5)
    assert row["clonotype_entropy"] == pytest.approx(
        -(2 / 3) * np.log(2 / 3) - (1 / 3) * np.log(1 / 3)
    )
    assert row["dominant_clonotype_fraction"] == pytest.approx(2 / 3)
    assert row["rfu_threshold_pass_rate"] == pytest.approx(2 / 3)
    assert row["cell_abundance"] == 3
    assert row["convergence_richness"] == 2
    assert row["mean_sequence_multiplicity"] == pytest.approx(3 / 2)
    assert row["normalized_convergence"] == pytest.approx(2 / 4)
    assert row["dominant_sequence_fraction"] == pytest.approx(2 / 3)
    assert row["threshold_pass_rate"] == pytest.approx(2 / 3)
    assert row["donor_count"] == 2
    assert row["donor_prevalence"] == pytest.approx(2 / 3)
    assert row["sample_count"] == 2
    assert row["sample_prevalence"] == pytest.approx(2 / 3)


def test_rfu_metrics_unique_sequence_weighting_is_explicit() -> None:
    result = rfu_metrics(_metrics_frame(), groupby=["phenotype"], weighting="unique_sequence")
    row = result[(result["phenotype"] == "P1") & (result["rfu_label"] == "RFU1")].iloc[0]

    assert row["weighting"] == "unique_sequence"
    assert row["weighted_abundance"] == pytest.approx(2 / 4)
    assert row["rfu_threshold_pass_rate"] == pytest.approx(0.5)


def test_rfu_metrics_prevalence_denominator_includes_unassigned_metadata() -> None:
    frame = pd.concat(
        [
            _metrics_frame(),
            pd.DataFrame(
                {
                    "cell_id": ["c7"],
                    "phenotype": ["P1"],
                    "cdr3aa": [pd.NA],
                    "rfu_label": [pd.NA],
                    "pass_thr": [pd.NA],
                    "donor": ["d5"],
                    "sample": ["s5"],
                }
            ),
        ],
        ignore_index=True,
    )
    result = rfu_metrics(
        frame,
        groupby="phenotype",
        weighting="cell",
        donor_col="donor",
        sample_col="sample",
    )
    row = result[(result["phenotype"] == "P1") & (result["rfu_label"] == "RFU1")].iloc[0]

    assert row["group_donor_count"] == 4
    assert row["donor_prevalence"] == pytest.approx(2 / 4)
    assert row["group_sample_count"] == 4
    assert row["sample_prevalence"] == pytest.approx(2 / 4)


@pytest.mark.parametrize("weighting", [None, "", "automatic", "cells"])
def test_rfu_metrics_rejects_implicit_or_unknown_weighting(weighting: object) -> None:
    with pytest.raises(ValueError, match="explicitly 'cell' or 'unique_sequence'"):
        rfu_metrics(_metrics_frame(), groupby="phenotype", weighting=weighting)  # type: ignore[arg-type]


def test_rfu_metrics_requires_explicit_grouping() -> None:
    with pytest.raises(ValueError, match="at least one explicit phenotype"):
        rfu_metrics(_metrics_frame(), groupby=[], weighting="cell")
