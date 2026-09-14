import numpy as np
import pandas as pd
import pytest

pytest.importorskip("sklearn")

from examples.regulatory_rfu_state_validation import (
    audit_predictions,
    paired_deltas,
    split_donor,
    training_support,
    weights,
)


def test_donor_split_purges_public_clones_and_all_timepoints():
    frame = pd.DataFrame(
        {"patient": ["A", "A", "B", "B", "C"], "clone_key": ["public", "a", "public", "b", "c"]}
    )
    train, test = split_donor(frame, "A")
    assert set(train.clone_key) == {"b", "c"}
    assert set(test.clone_key) == {"public", "a"}
    assert not set(train.patient) & set(test.patient)
    ordinary_train, ordinary_test = split_donor(frame, "A", purge_shared=False)
    assert "public" in set(ordinary_train.clone_key)
    pd.testing.assert_frame_equal(test, ordinary_test)


def test_paired_log_loss_delta_sign_and_prediction_order():
    rows, metrics = [], []
    for name, true_probability in (("TRBV_TRBJ_length", 0.5), ("TRBV_TRBJ_length_RFU", 0.8)):
        rows.append(
            {
                "cell_id": "cell",
                "held_out_donor": "A",
                "weighting": "clone",
                "model": name,
                "cluster": "z",
                "evaluation_weight": 1.0,
                "prob:z": true_probability,
                "prob:a": 1 - true_probability,
            }
        )
        metrics.append(
            {
                "held_out_donor": "A",
                "weighting": "clone",
                "model": name,
                "log_loss": -np.log(true_probability),
            }
        )
    metrics = pd.DataFrame(metrics)
    predictions = pd.DataFrame(rows)
    audit_predictions(predictions, metrics)
    assert paired_deltas(metrics).delta_log_loss.iloc[0] == pytest.approx(np.log(0.5 / 0.8))
    metrics.loc[metrics.model.str.endswith("_RFU"), "log_loss"] = 2.0
    assert paired_deltas(metrics).delta_log_loss.iloc[0] > 0
    with pytest.raises(ValueError, match="disagree"):
        audit_predictions(predictions, metrics)


def test_prediction_audit_rejects_different_evaluation_cells():
    rows = pd.DataFrame(
        [
            {
                "cell_id": c,
                "cluster": "z",
                "evaluation_weight": 1.0,
                "held_out_donor": "A",
                "weighting": "clone",
                "model": m,
            }
            for c, m in [("c1", "TRBV_TRBJ_length"), ("c2", "TRBV_TRBJ_length_RFU")]
        ]
    )
    with pytest.raises(AssertionError):
        audit_predictions(rows, pd.DataFrame())


def test_expansion_does_not_change_total_clone_weight():
    frame = pd.DataFrame(
        {"patient": ["A"] * 4 + ["B"], "clone_key": ["expanded"] * 3 + ["single", "other"]}
    )
    w = weights(frame, "clone", balance_donors=False)
    assert np.isclose(w[:3].sum(), w[3])
    balanced = weights(frame, "clone", balance_donors=True)
    assert np.isclose(balanced[:4].sum(), balanced[4])


def test_training_support_counts_clones_not_cells_or_test_donors():
    train = pd.DataFrame(
        {
            "patient": ["A"] * 100 + ["B", "C", "A"],
            "clone_key": ["a"] * 100 + ["b", "c", "x"],
            "rfu": ["supported"] * 102 + ["rare"],
        }
    )
    assert training_support(train, "rfu", 3) == {"supported"}
    assert training_support(train, "rfu", 4) == set()
