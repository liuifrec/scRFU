from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from scrfu.tl import (
    annotate_vdjdb,
    compare_antigen_groupings,
    global_antigen_coherence,
    load_vdjdb_reference,
    rfu_antigen_abundance,
    rfu_antigen_coherence,
    rfu_antigen_permutation_test,
    summarize_antigen_context,
)


def _results() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "input_row_id": ["r1", "r2", "r3", "r4", "r5", "r6"],
            "unique_sequence_id": ["s1", "s1", "s2", "s3", "s4", "s5"],
            "cell_id": ["c1", "c2", "c3", "c4", "c5", "c6"],
            "cdr3aa": ["CASSAAA", "CASSAAA", "CASSAAD", "CASSCCC", "CASSCCD", "CASSGGG"],
            "v_call": ["TRBV1", "TRBV1", "TRBV1", "TRBV2", "TRBV2", "TRBV3"],
            "chain": ["TRB"] * 6,
            "rfu_label": ["R1", "R1", "R1", "R2", "R2", "R3"],
            "pass_thr": [True, True, True, True, False, True],
        }
    )


def _evidence() -> pd.DataFrame:
    reference = load_vdjdb_reference(
        pd.DataFrame(
            {
                "cdr3aa": ["CASSAAA", "CASSAAD", "CASSCCC", "CASSCCD", "CASSCCD"],
                "v_call": ["TRBV1", "TRBV1", "TRBV2", "TRBV2", "TRBV2"],
                "chain": ["TRB"] * 5,
                "epitope": ["A", "A", "A", "A", "B"],
                "score": [3, 2, 2, 1, 2],
            }
        ),
        release_label="synthetic",
    )
    return annotate_vdjdb(_results(), reference)


def test_rfu_coherence_hand_calculated_sequence_weighting() -> None:
    result = rfu_antigen_coherence(_results(), _evidence(), min_matched_sequences=2)
    r1 = result.set_index("rfu_label").loc["R1"]
    assert r1["total_rfu_sequences"] == 2
    assert r1["total_rfu_cells"] == 3
    assert r1["vdjdb_matched_sequences"] == 2
    assert r1["vdjdb_matched_cells"] == 3
    assert r1["sequence_match_rate"] == 1
    assert r1["cell_match_rate"] == 1
    assert r1["antigen_purity"] == 1
    assert r1["antigen_entropy"] == 0
    assert r1["normalized_antigen_entropy"] == 0
    r2 = result.set_index("rfu_label").loc["R2"]
    assert r2["antigen_purity"] == pytest.approx(0.75)
    assert r2["antigen_entropy"] == pytest.approx(-(0.75 * math.log(0.75) + 0.25 * math.log(0.25)))
    assert 0 <= r2["normalized_antigen_entropy"] <= 1
    assert r2["ambiguous_sequence_fraction"] == pytest.approx(0.5)
    assert not result.set_index("rfu_label").loc["R3", "eligible_for_coherence"]


def test_ambiguity_policies_and_cell_weighting() -> None:
    fractional = rfu_antigen_abundance(_results(), _evidence(), ambiguity_policy="fractional")
    r2 = fractional.loc[fractional["rfu_label"].eq("R2")].set_index("epitope")
    assert r2.loc["A", "antigen_abundance"] == pytest.approx(1.5)
    assert r2.loc["B", "antigen_abundance"] == pytest.approx(0.5)
    excluded = rfu_antigen_abundance(_results(), _evidence(), ambiguity_policy="exclude_ambiguous")
    assert excluded.loc[excluded["rfu_label"].eq("R2"), "epitope"].tolist() == ["A"]
    cells = rfu_antigen_coherence(
        _results(), _evidence(), weighting="cell", min_matched_sequences=1
    )
    assert cells.set_index("rfu_label").loc["R1", "antigen_purity"] == 1
    cell_abundance = rfu_antigen_abundance(_results(), _evidence(), weighting="cell")
    assert cell_abundance.set_index(["rfu_label", "epitope"]).loc[
        ("R1", "A"), "antigen_abundance"
    ] == pytest.approx(3)


def test_threshold_policy_changes_selected_groups() -> None:
    nearest = rfu_antigen_coherence(_results(), _evidence(), min_matched_sequences=1)
    threshold = rfu_antigen_coherence(
        _results(), _evidence(), assignment_policy="threshold_pass", min_matched_sequences=1
    )
    assert nearest.set_index("rfu_label").loc["R2", "total_rfu_sequences"] == 2
    assert threshold.set_index("rfu_label").loc["R2", "total_rfu_sequences"] == 1


def test_global_coherence_and_bounds() -> None:
    result = global_antigen_coherence(_results(), _evidence())
    assert 0 <= result["weighted_mean_rfu_antigen_purity"] <= 1
    assert 0 <= result["same_antigen_pair_fraction_within_rfus"] <= 1
    assert 0 <= result["between_rfu_antigen_sharing"] <= 1
    assert result["mutual_information"] >= 0
    assert 0 <= result["normalized_mutual_information"] <= 1


def test_permutation_is_deterministic_preserves_sizes_and_does_not_mutate() -> None:
    original = _results()
    before = original.copy(deep=True)
    first = rfu_antigen_permutation_test(original, _evidence(), n_permutations=20, random_state=7)
    second = rfu_antigen_permutation_test(original, _evidence(), n_permutations=20, random_state=7)
    np.testing.assert_array_equal(first.permutation_values, second.permutation_values)
    assert first.parameters["group_sizes_preserved"]
    assert first.empirical_upper_tail_probability == pytest.approx(
        (1 + np.sum(first.permutation_values >= first.observed)) / 21
    )
    pd.testing.assert_frame_equal(original, before)


@pytest.mark.parametrize("stratify_by", ["cdr3_length", "TRBV", ["cdr3_length", "trbv"]])
def test_stratified_permutation(stratify_by) -> None:
    result = rfu_antigen_permutation_test(
        _results(), _evidence(), n_permutations=5, random_state=3, stratify_by=stratify_by
    )
    assert len(result.permutation_values) == 5


def test_permutation_validation_and_zero_variance() -> None:
    with pytest.raises(ValueError, match="positive integer"):
        rfu_antigen_permutation_test(_results(), _evidence(), n_permutations=0)
    one = _results().loc[_results()["unique_sequence_id"].eq("s1")]
    one_evidence = _evidence().loc[_evidence()["unique_sequence_id"].eq("s1")]
    with pytest.raises(ValueError, match="At least two"):
        rfu_antigen_permutation_test(one, one_evidence, n_permutations=2)
    result = rfu_antigen_permutation_test(
        _results(),
        _evidence(),
        n_permutations=5,
        random_state=0,
        stratify_by=["cdr3_length", "trbv"],
    )
    if result.null_std == 0:
        assert math.isnan(result.z_score)


def test_grouping_baselines_are_tidy_and_deterministic() -> None:
    first = compare_antigen_groupings(
        _results(),
        _evidence(),
        groupings=("rfu", "trbv", "cdr3_length", "trbv_cdr3_length", "size_matched_random"),
        random_state=4,
    )
    second = compare_antigen_groupings(
        _results(),
        _evidence(),
        groupings=("rfu", "trbv", "cdr3_length", "trbv_cdr3_length", "size_matched_random"),
        random_state=4,
    )
    pd.testing.assert_frame_equal(first, second)
    assert set(first["metric"]) == {"purity", "entropy", "same_antigen_pair_fraction"}


def _baseline_case(*, rfu_coherent: bool) -> tuple[pd.DataFrame, pd.DataFrame]:
    results = pd.DataFrame(
        {
            "input_row_id": ["a", "b", "c", "d"],
            "unique_sequence_id": ["a", "b", "c", "d"],
            "cdr3aa": ["CASSAAA", "CASSAAD", "CASSCCC", "CASSCCD"],
            "v_call": ["TRBV1", "TRBV2", "TRBV1", "TRBV2"]
            if rfu_coherent
            else ["TRBV1", "TRBV2", "TRBV1", "TRBV2"],
            "chain": ["TRB"] * 4,
            "rfu_label": ["R1", "R1", "R2", "R2"],
            "pass_thr": [True] * 4,
        }
    )
    labels = ["A", "A", "B", "B"] if rfu_coherent else ["A", "B", "A", "B"]
    reference = load_vdjdb_reference(
        results[["cdr3aa", "v_call", "chain"]].assign(epitope=labels),
        release_label="baseline",
    )
    return results, annotate_vdjdb(results, reference)


def test_grouping_baselines_can_favor_rfu_or_trbv() -> None:
    rfu_results, rfu_evidence = _baseline_case(rfu_coherent=True)
    rfu_case = compare_antigen_groupings(
        rfu_results, rfu_evidence, groupings=("rfu", "trbv"), metrics=("purity",)
    ).set_index("grouping_method")
    assert rfu_case.loc["rfu", "value"] > rfu_case.loc["trbv", "value"]
    global_rfu = global_antigen_coherence(rfu_results, rfu_evidence)
    assert global_rfu["same_antigen_pair_fraction_within_rfus"] == 1
    assert global_rfu["mutual_information"] == pytest.approx(math.log(2))
    assert global_rfu["normalized_mutual_information"] == pytest.approx(1)

    trbv_results, trbv_evidence = _baseline_case(rfu_coherent=False)
    trbv_case = compare_antigen_groupings(
        trbv_results, trbv_evidence, groupings=("rfu", "trbv"), metrics=("purity",)
    ).set_index("grouping_method")
    assert trbv_case.loc["trbv", "value"] > trbv_case.loc["rfu", "value"]
    mixed = rfu_antigen_coherence(trbv_results, trbv_evidence, min_matched_sequences=2).set_index(
        "rfu_label"
    )
    assert mixed.loc["R1", "antigen_purity"] == pytest.approx(0.5)
    assert mixed.loc["R1", "normalized_antigen_entropy"] == pytest.approx(1.0)


def test_one_matched_sequence_is_not_coherence_eligible() -> None:
    results = _results().loc[_results()["rfu_label"].eq("R1")].copy()
    evidence = _evidence().loc[_evidence()["unique_sequence_id"].eq("s1")]
    result = rfu_antigen_coherence(results, evidence, min_matched_sequences=2).iloc[0]
    assert result["sequence_match_rate"] == pytest.approx(0.5)
    assert not result["eligible_for_coherence"]


def test_metadata_context_alignment_prevalence_and_recurrence() -> None:
    metadata = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4", "c5", "c6"],
            "sample": ["a", "b", "a", "a", "b", "b"],
            "donor": ["d1", "d2", "d1", "d1", "d2", "d2"],
            "phenotype": ["T1", "T1", "T2", "T2", "T2", "T1"],
        }
    )
    result = summarize_antigen_context(
        _results(),
        _evidence(),
        metadata=metadata,
        sample_key="sample",
        donor_key="donor",
        phenotype_key="phenotype",
    )
    assert result.join_qc["input_row_count"] == 6
    assert result.join_qc["output_row_count"] == 6
    assert (
        result.sequence_prevalence.set_index("unique_sequence_id").loc[
            "s1", "sample_prevalence_count"
        ]
        == 2
    )
    assert not result.phenotype_abundance.empty
    duplicate = pd.concat([metadata, metadata.iloc[[0]]])
    with pytest.raises(ValueError, match="unique"):
        summarize_antigen_context(_results(), _evidence(), metadata=duplicate)
    missing = summarize_antigen_context(
        _results(),
        _evidence(),
        metadata=metadata.iloc[:-1],
        sample_key="sample",
    )
    assert missing.join_qc["unmatched_metadata_rows"] == 1
    conflicting = metadata.copy()
    conflicting.loc[conflicting["cell_id"].eq("c2"), "donor"] = "d1"
    with pytest.raises(ValueError, match="Conflicting"):
        summarize_antigen_context(
            _results(),
            _evidence(),
            metadata=conflicting,
            sample_key="sample",
            donor_key="donor",
        )
