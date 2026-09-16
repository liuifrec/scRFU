from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from scrfu.tl import multiscale_repertoire_change, permute_fixed_groups


def series(values):
    return pd.Series(values, dtype=float)


def test_changes_within_group_are_hidden_but_across_groups_retained():
    before = series({"a": 10, "b": 0, "c": 0})
    within = series({"a": 0, "b": 10, "c": 0})
    across = series({"a": 0, "b": 0, "c": 10})
    mapping = pd.Series({"a": "g1", "b": "g1", "c": "g2"})
    x = multiscale_repertoire_change(before, within, mapping).summary
    assert (x["d_clone"], x["d_group"], x["aggregation_cancellation"]) == (1, 0, 1)
    y = multiscale_repertoire_change(before, across, mapping).summary
    assert (y["d_clone"], y["d_group"], y["aggregation_cancellation"]) == (1, 1, 0)


def test_one_clone_expansion_and_scale_invariance():
    a, b = series({"a": 1, "b": 1}), series({"a": 9, "b": 1})
    mapping = pd.Series({"a": "g1", "b": "g2"})
    for left, right in [(a, b), (a / a.sum(), b / b.sum())]:
        x = multiscale_repertoire_change(left, right, mapping)
        assert x.summary["d_clone"] == pytest.approx(0.4)
        assert x.summary["d_group"] == pytest.approx(0.4)
        assert x.summary["aggregation_cancellation"] == 0


def test_changed_coverage_is_not_hidden_and_missing_visits_rejected():
    a, b = series({"mapped": 8, "unmapped": 2}), series({"mapped": 2, "unmapped": 8})
    mapping = pd.Series({"mapped": "g"})
    with pytest.raises(ValueError, match="conditioning"):
        multiscale_repertoire_change(a, b, mapping)
    x = multiscale_repertoire_change(a, b, mapping, unmapped="condition").summary
    assert x["d_clone"] == x["d_group"] == 0
    assert x["coverage_before"] == 0.8 and x["coverage_after"] == 0.2
    assert x["estimand"] == "conditional_on_mapped_receptors"
    with pytest.raises(TypeError, match="missing visits"):
        multiscale_repertoire_change(a, None, mapping)


def test_empty_is_unknown_identical_nonempty_is_zero():
    empty, a = series({}), series({"a": 2})
    mapping = pd.Series({"a": "g"})
    for before, after in [(empty, empty), (a, empty)]:
        x = multiscale_repertoire_change(before, after, mapping).summary
        assert x["status"] == "empty_mapped_sample"
        assert np.isnan(x["d_clone"]) and np.isnan(x["d_group"])
    assert multiscale_repertoire_change(a, a, mapping).summary["d_clone"] == 0


@pytest.mark.parametrize("bad", [-1, np.nan, np.inf])
def test_invalid_mass(bad):
    with pytest.raises(ValueError, match="finite and nonnegative"):
        multiscale_repertoire_change(series({"a": bad}), series({"a": 1}), pd.Series({"a": "g"}))


def test_mapping_cannot_change_between_visits_or_duplicate_receptors():
    a = series({"a": 1})
    with pytest.raises(ValueError, match="unique"):
        multiscale_repertoire_change(a, a, pd.Series(["g1", "g2"], index=["a", "a"]))
    with pytest.raises(TypeError, match="Series"):
        multiscale_repertoire_change(a, a, pd.DataFrame({"before": ["g1"], "after": ["g2"]}))


def test_sampling_only_does_not_imply_latent_change_and_contraction_holds():
    rng = np.random.default_rng(7)
    names = [f"c{i}" for i in range(40)]
    mapping = pd.Series([f"g{i // 4}" for i in range(40)], index=names)
    latent = np.repeat(1 / 40, 40)
    positive = 0
    for _ in range(30):
        a = pd.Series(rng.multinomial(100, latent), index=names)
        b = pd.Series(rng.multinomial(100, latent), index=names)
        x = multiscale_repertoire_change(a, b, mapping).summary
        assert 0 <= x["d_group"] <= x["d_clone"] + 1e-12 <= 1 + 1e-12
        positive += x["d_clone"] > 0
    assert positive == 30  # Simulated observation noise despite identical latent frequencies.


def test_fixed_random_groups_preserve_sizes_features_and_input_order_invariance():
    names = [f"c{i:02}" for i in range(20)]
    mapping = pd.Series(["g1", "g2"] * 10, index=names)
    strata = pd.DataFrame({"v": ["V1"] * 10 + ["V2"] * 10}, index=names)
    perm = permute_fixed_groups(mapping, strata=strata, random_state=10)
    reverse = permute_fixed_groups(mapping.iloc[::-1], strata=strata.iloc[::-1], random_state=10)
    pd.testing.assert_series_equal(perm.sort_index(), reverse.sort_index())
    for feature in ("V1", "V2"):
        idx = strata.index[strata.v.eq(feature)]
        assert mapping.loc[idx].value_counts().to_dict() == perm.loc[idx].value_counts().to_dict()
    assert perm.ne(mapping).any()


def test_unexchangeable_random_strata_stay_fixed_and_missing_strata_fail():
    mapping = pd.Series({"a": "g1", "b": "g2"})
    strata = pd.DataFrame({"v": ["V1", "V2"]}, index=mapping.index)
    pd.testing.assert_series_equal(mapping, permute_fixed_groups(mapping, strata=strata))
    with pytest.raises(ValueError, match="matching strata"):
        permute_fixed_groups(mapping, strata=strata.iloc[:1])
