from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import scrfu


def _rows() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "input_row_id": [f"r{i}" for i in range(8)],
            "cell_id": [f"c{i}" for i in range(8)],
            "cdr3aa": [f"CASS{i}F" for i in range(8)],
            "rfu_label": ["RFU1", "RFU2", "RFU1", "RFU3", "RFU2", "RFU3", "RFU1", "RFU3"],
            "pass_thr": [True, True, True, True, True, False, True, True],
            "sample": ["d1_t0", "d1_t0", "d1_t1", "d1_t1", "d2_t0", "d2_t0", "d2_t1", "d2_t1"],
            "donor": ["d1", "d1", "d1", "d1", "d2", "d2", "d2", "d2"],
            "time": [0, 0, 1, 1, 0, 0, 1, 1],
            "compartment": ["A", "A", "A", "A", "A", "A", "A", "A"],
        }
    )


def test_validate_longitudinal_design_balanced_unbalanced_and_invalid() -> None:
    balanced = scrfu.tl.validate_longitudinal_design(
        _rows(), sample_key="sample", donor_key="donor", time_key="time"
    )
    assert balanced.ordered_donors == ("d1", "d2")
    assert balanced.ordered_timepoints == (0, 1)
    assert balanced.qc_table["repeated_donor"].all()

    unbalanced_rows = _rows().loc[lambda x: x["sample"].ne("d2_t1")]
    unbalanced = scrfu.tl.validate_longitudinal_design(
        unbalanced_rows, sample_key="sample", donor_key="donor", time_key="time"
    )
    d2 = unbalanced.qc_table.set_index("donor").loc["d2"]
    assert d2["singleton_donor"]
    assert d2["missing_timepoint_count"] == 1

    conflicting = _rows()
    conflicting.loc[1, "donor"] = "wrong"
    with pytest.raises(ValueError, match="conflicting"):
        scrfu.tl.validate_longitudinal_design(
            conflicting, sample_key="sample", donor_key="donor", time_key="time"
        )
    ambiguous = _rows().assign(time=lambda x: x["time"].map({0: "baseline", 1: "later"}))
    with pytest.raises(ValueError, match="Ambiguous time labels"):
        scrfu.tl.validate_longitudinal_design(
            ambiguous, sample_key="sample", donor_key="donor", time_key="time"
        )


@pytest.mark.parametrize("normalize", ["count", "proportion", "counts_per_1000", "clr"])
def test_longitudinal_matrix_normalization_and_missingness(normalize: str) -> None:
    rows = _rows().loc[lambda x: x["sample"].ne("d2_t1")]
    result = scrfu.tl.rfu_longitudinal_matrix(
        rows,
        sample_key="sample",
        donor_key="donor",
        time_key="time",
        normalize=normalize,
    )
    assert result.sample_matrix.index.tolist() == ["d1_t0", "d1_t1", "d2_t0"]
    assert result.parameters["no_imputation"] is True
    assert result.missingness_mask.loc["d2"].to_numpy().sum() == 1
    if normalize == "proportion":
        np.testing.assert_allclose(result.sample_matrix.sum(axis=1), 1)
    elif normalize == "counts_per_1000":
        np.testing.assert_allclose(result.sample_matrix.sum(axis=1), 1000)
    elif normalize == "clr":
        np.testing.assert_allclose(result.sample_matrix.mean(axis=1), 0, atol=1e-12)


@pytest.mark.parametrize(
    ("metric", "direction"),
    [
        ("cosine", "similarity"),
        ("jaccard", "similarity"),
        ("weighted_jaccard", "similarity"),
        ("bray_curtis_dissimilarity", "distance"),
        ("jensen_shannon_distance", "distance"),
    ],
)
def test_longitudinal_similarity_tidy_pairs(metric: str, direction: str) -> None:
    result = scrfu.tl.rfu_longitudinal_matrix(
        _rows(), sample_key="sample", donor_key="donor", time_key="time"
    )
    pairs = scrfu.tl.longitudinal_similarity(result, metric=metric)
    assert len(pairs) == 6
    assert pairs["status"].eq("valid").all()
    assert pairs["direction"].eq(direction).all()
    assert set(pairs["time_interval"]) == {0.0, 1.0}
    summary = scrfu.tl.summarize_longitudinal_similarity(pairs)
    assert set(summary["donor_relation"]) == {"within_donor", "between_donor"}


def test_donor_retrieval_recoverable_and_absent() -> None:
    matrix = pd.DataFrame(
        [[5, 0], [4, 0], [0, 5], [0, 4]],
        index=["d1_t0", "d1_t1", "d2_t0", "d2_t1"],
        columns=["RFU1", "RFU2"],
    )
    metadata = pd.DataFrame(
        {"donor": ["d1", "d1", "d2", "d2"], "time": [0, 1, 0, 1]},
        index=matrix.index,
    )
    result = scrfu.tl.donor_retrieval(matrix, metadata=metadata, top_k=1)
    assert result["top_1_donor_match"].all()
    assert result.attrs["mean_reciprocal_rank"] == 1

    singleton = scrfu.tl.donor_retrieval(
        matrix.iloc[[0, 2]], metadata=metadata.iloc[[0, 2]], top_k=1
    )
    assert not singleton["top_k_donor_match"].any()
    assert singleton["correct_donor_rank"].isna().all()


def test_dynamics_compartment_and_cluster_resampling() -> None:
    rows = _rows()
    longitudinal = scrfu.tl.rfu_longitudinal_matrix(
        rows, sample_key="sample", donor_key="donor", time_key="time"
    )
    dynamics = scrfu.tl.rfu_longitudinal_dynamics(
        longitudinal, expansion_fold_change=2, contraction_fold_change=2
    )
    labels = dynamics.classifications.set_index(["donor", "rfu_label"])["classification"]
    assert labels.loc[("d1", "RFU3")] == "appearing"
    assert labels.loc[("d2", "RFU2")] == "disappearing"

    compartment_rows = pd.concat(
        [
            rows.assign(compartment="A"),
            rows.assign(compartment="B", sample=lambda x: x["sample"] + "_b"),
        ],
        ignore_index=True,
    )
    compartment = scrfu.tl.rfu_longitudinal_matrix(
        compartment_rows,
        sample_key="sample",
        donor_key="donor",
        time_key="time",
        compartment_key="compartment",
    )
    comparison = scrfu.tl.longitudinal_compartment_comparison(compartment)
    assert len(comparison.pairwise_comparison) == 4
    assert set(comparison.rfu_status["compartment_status"]) == {"shared", "absent"}

    donor_data = pd.DataFrame(
        {"donor": ["d1", "d1", "d2", "d2"], "label": ["A", "A", "B", "B"], "x": [1, 2, 3, 4]}
    )

    def statistic(frame: pd.DataFrame) -> float:
        return float(frame["x"].mean())

    boot_a = scrfu.tl.bootstrap_longitudinal_statistic(
        donor_data, statistic, donor_key="donor", n_resamples=10, random_state=4
    )
    boot_b = scrfu.tl.bootstrap_longitudinal_statistic(
        donor_data, statistic, donor_key="donor", n_resamples=10, random_state=4
    )
    pd.testing.assert_frame_equal(boot_a.values, boot_b.values)
    permuted = scrfu.tl.permute_longitudinal_labels(
        donor_data,
        lambda frame: float(frame.loc[frame["label"].eq("A"), "x"].mean()),
        donor_key="donor",
        label_key="label",
        exact=True,
    )
    assert permuted.parameters["actual_permutations"] == 2
    assert 0 <= permuted.summary["empirical_upper_tail_probability"] <= 1
