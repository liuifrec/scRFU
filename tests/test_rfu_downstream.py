from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from scrfu.tl import rfu_overlap, rfu_phenotype_coupling, rfu_pseudobulk


def _assignments() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cell_id": ["a1", "a2", "a3", "b1", "b2", "c1", "c2"],
            "sample": ["s2", "s2", "s2", "s1", "s1", "s3", "s3"],
            "phenotype": ["T", "T", "T", "B", "B", "T", "T"],
            "cdr3aa": ["A", "A", "B", "A", "C", "D", "E"],
            "rfu_label": ["RFU2", "RFU2", "RFU10", "RFU2", "RFU10", "RFU2", pd.NA],
            "pass_thr": [True, False, True, True, True, True, False],
        }
    )


def test_pseudobulk_counts_normalizations_and_order() -> None:
    counts = rfu_pseudobulk(
        _assignments(), sample_key="sample", phenotype_keys="phenotype", weighting="cell"
    )
    assert counts.matrix.index.tolist() == ["s1", "s2", "s3"]
    assert counts.matrix.columns.tolist() == ["RFU2", "RFU10"]
    assert counts.matrix.loc["s2"].tolist() == [2, 1]
    assert counts.sample_metadata.loc["s1", "phenotype"] == "B"

    proportions = rfu_pseudobulk(
        _assignments(), sample_key="sample", weighting="cell", normalize="proportion"
    ).matrix
    per_1000 = rfu_pseudobulk(
        _assignments(), sample_key="sample", weighting="cell", normalize="counts_per_1000"
    ).matrix
    assert proportions.loc["s2", "RFU2"] == pytest.approx(2 / 3)
    assert per_1000.loc["s2", "RFU2"] == pytest.approx(2000 / 3)

    clr = rfu_pseudobulk(
        _assignments(), sample_key="sample", weighting="cell", normalize="clr", pseudocount=1
    ).matrix
    expected = math.log(3) - (math.log(3) + math.log(2)) / 2
    assert clr.loc["s2", "RFU2"] == pytest.approx(expected)
    assert clr.sum(axis=1).to_numpy() == pytest.approx(np.zeros(3))


def test_pseudobulk_policy_sequence_weighting_prevalence_and_empty_sample() -> None:
    threshold = rfu_pseudobulk(
        _assignments(), sample_key="sample", assignment_policy="threshold_pass", weighting="cell"
    )
    assert threshold.matrix.loc["s2", "RFU2"] == 1
    sequence = rfu_pseudobulk(_assignments(), sample_key="sample", weighting="unique_sequence")
    assert sequence.matrix.loc["s2", "RFU2"] == 1
    filtered = rfu_pseudobulk(
        _assignments(), sample_key="sample", weighting="cell", min_prevalence=3
    )
    assert filtered.matrix.columns.tolist() == ["RFU2"]
    assert filtered.matrix.loc["s3", "RFU2"] == 1

    categorical = _assignments()
    categorical["rfu_label"] = pd.Categorical(
        categorical["rfu_label"], categories=["RFU2", "RFU10", "RFU99"]
    )
    retained = rfu_pseudobulk(
        categorical, sample_key="sample", retain_zero_rfus=True, weighting="cell"
    )
    assert retained.matrix["RFU99"].eq(0).all()


def test_pseudobulk_rejects_sample_metadata_conflicts() -> None:
    frame = _assignments()
    frame.loc[1, "phenotype"] = "B"
    with pytest.raises(ValueError, match="conflicting"):
        rfu_pseudobulk(frame, sample_key="sample", phenotype_keys="phenotype")


@pytest.mark.parametrize(
    ("metric", "expected", "direction"),
    [
        ("jaccard", 1 / 2, "similarity"),
        ("sorensen_dice", 2 / 3, "similarity"),
        ("overlap_coefficient", 1.0, "similarity"),
        ("cosine", 1 / math.sqrt(5), "similarity"),
        ("bray_curtis_similarity", 0.5, "similarity"),
        ("bray_curtis_dissimilarity", 0.5, "distance"),
        ("weighted_jaccard", 1 / 3, "similarity"),
    ],
)
def test_overlap_hand_calculated_metrics(metric: str, expected: float, direction: str) -> None:
    matrix = pd.DataFrame([[1, 0], [1, 2]], index=["a", "b"], columns=["x", "y"])
    result = rfu_overlap(matrix, metric=metric)
    assert result.matrix.loc["a", "b"] == pytest.approx(expected)
    assert result.direction == direction


def test_overlap_jensen_shannon_and_zero_vectors() -> None:
    identical = pd.DataFrame([[1, 1], [2, 2]], index=["a", "b"])
    result = rfu_overlap(identical, metric="jensen_shannon_distance")
    assert result.matrix.loc["a", "b"] == pytest.approx(0.0)
    zeros = pd.DataFrame([[0, 0], [1, 0]], index=["zero", "nonzero"])
    result = rfu_overlap(zeros, metric="jaccard")
    assert math.isnan(result.matrix.loc["zero", "zero"])
    assert result.matrix.loc["zero", "nonzero"] == 0


def _coupling() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cell_id": ["a", "b", "c", "d", "e", "f"],
            "sample": ["s1", "s2", "s3", "s1", "s2", "s3"],
            "phenotype": ["P1", "P1", "P1", "P2", "P2", "P2"],
            "cdr3aa": ["A", "B", "C", "D", "E", "F"],
            "rfu_label": ["one", "one", "one", "even", "even", "even"],
            "pass_thr": [True, False, True, True, True, False],
        }
    )


def test_phenotype_coupling_entropy_specificity_and_recurrence() -> None:
    result = rfu_phenotype_coupling(
        _coupling(), phenotype_key="phenotype", sample_key="sample", weighting="cell"
    )
    one = result.loc[result["rfu_label"].eq("one")].iloc[0]
    assert one["phenotype_entropy"] == 0
    assert one["normalized_phenotype_entropy"] == 0
    assert one["phenotype_specificity"] == 1
    assert one["dominant_phenotype"] == "P1"
    assert one["sample_recurrence_count"] == 3
    even = result.loc[result["rfu_label"].eq("even")].iloc[0]
    assert even["phenotype_entropy"] == 0
    assert even["phenotype_specificity"] == 1


def test_phenotype_coupling_even_distribution_and_threshold() -> None:
    frame = _coupling()
    frame.loc[2, "rfu_label"] = "shared"
    frame.loc[5, "rfu_label"] = "shared"
    result = rfu_phenotype_coupling(frame, phenotype_key="phenotype", weighting="cell")
    shared = result.loc[result["rfu_label"].eq("shared")].iloc[0]
    assert shared["normalized_phenotype_entropy"] == pytest.approx(1)
    assert shared["phenotype_specificity"] == pytest.approx(0)
    threshold = rfu_phenotype_coupling(
        frame,
        phenotype_key="phenotype",
        assignment_policy="threshold_pass",
        weighting="unique_sequence",
    )
    assert "shared" in threshold["rfu_label"].values
