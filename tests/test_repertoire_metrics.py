from __future__ import annotations

import math

import pandas as pd
import pytest

from scrfu.tl import repertoire_metrics


def _repertoire() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample": ["s1"] * 5 + ["s2"],
            "cell_id": ["c1", "c2", "c3", "c4", "c5", "c6"],
            "cdr3aa": ["AAA", "AAA", "BBB", "CCC", "CCC", "DDD"],
            "clonotype_id": ["cl1", "cl1", "cl2", "cl3", "cl3", "cl4"],
            "v_call": ["TRBV1", "TRBV1", "TRBV2", "TRBV2", "TRBV2", "TRBV3"],
            "productive": [True, True, True, False, True, True],
            "chain": ["TRB"] * 6,
        }
    )


def test_repertoire_metrics_hand_calculated_cell_weighting() -> None:
    result = repertoire_metrics(_repertoire(), groupby="sample", weighting="cell")
    row = result.loc[result["sample"].eq("s1")].iloc[0]
    probabilities = [2 / 5, 1 / 5, 2 / 5]
    assert row["receptor_bearing_row_count"] == 5
    assert row["unique_cell_count"] == 5
    assert row["unique_cdr3_richness"] == 3
    assert row["clonotype_richness"] == 3
    assert row["shannon_entropy"] == pytest.approx(-sum(p * math.log(p) for p in probabilities))
    assert row["simpson_diversity"] == pytest.approx(1 - sum(p * p for p in probabilities))
    assert row["inverse_simpson_diversity"] == pytest.approx(1 / sum(p * p for p in probabilities))
    assert row["dominant_cdr3_fraction"] == pytest.approx(2 / 5)
    assert row["dominant_clonotype_fraction"] == pytest.approx(2 / 5)
    assert row["productive_fraction"] == pytest.approx(4 / 5)
    assert row["trbv_richness"] == 2
    assert row["top_trbv_fraction"] == pytest.approx(3 / 5)


def test_repertoire_metrics_sequence_fallback_and_weighting() -> None:
    table = _repertoire().drop(columns="clonotype_id")
    result = repertoire_metrics(table, groupby="sample", weighting="unique_sequence")
    row = result.loc[result["sample"].eq("s1")].iloc[0]
    assert row["clonotype_fallback"]
    assert pd.isna(row["clonotype_richness"])
    assert row["diversity_identity"] == "cdr3aa_sequence"
    assert row["shannon_entropy"] == pytest.approx(math.log(3))
    assert row["dominant_clonotype_fraction"] is pd.NA


def test_repertoire_metrics_empty_and_explicit_semantics() -> None:
    empty = _repertoire().iloc[:0]
    assert repertoire_metrics(empty, groupby="sample", weighting="cell").empty
    with pytest.raises(ValueError, match="weighting"):
        repertoire_metrics(_repertoire(), groupby="sample", weighting="automatic")
    with pytest.raises(ValueError, match="explicit group"):
        repertoire_metrics(_repertoire(), groupby=[], weighting="cell")
