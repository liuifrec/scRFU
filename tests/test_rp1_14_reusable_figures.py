"""Figure contracts: prevent denominator, pairing and frozen-input mistakes."""

import json

import numpy as np
import pandas as pd
import pytest

from manuscript.scripts.rp1_14_reusable_figures import (
    build_tables,
    column_note,
    primary_rows,
    replicate_summary,
    safe_child,
    sha256,
    summarize_donors,
    verify_upstream,
)


def test_figure_runner_cannot_retune_endpoint_definitions():
    with pytest.raises(ValueError, match="selection is frozen"):
        build_tables({}, {"selection": {"threshold": 0.7}})


def test_figure_selection_keeps_frozen_endpoints_policy_and_weights():
    base = dict(
        donor="synthetic",
        compartment="CD4",
        primary_pair=True,
        policy="threshold",
        weighting="reads",
        grouping="RFU",
        status="valid",
        visit_before=1,
        visit_after=3,
        d_group=0.3,
    )
    frame = pd.DataFrame(
        [
            base,
            {**base, "primary_pair": False, "visit_after": 2},
            {**base, "policy": "nearest"},
            {**base, "weighting": "unique_clone"},
        ]
    )
    assert primary_rows(frame).index.tolist() == [0]
    assert primary_rows(frame, weighting="unique_clone").index.tolist() == [3]
    with pytest.raises(ValueError, match="Duplicate"):
        primary_rows(pd.concat([frame, frame.iloc[:1]]))
    with pytest.raises(ValueError, match="endpoints"):
        primary_rows(pd.DataFrame([{**base, "visit_after": 2}]))


def synthetic_replicates():
    return pd.DataFrame(
        [
            dict(
                donor=donor,
                compartment=comp,
                replicate=i,
                primary_pair=True,
                status="valid",
                d_clone=0.9,
                d_group=value,
                aggregation_cancellation=0.9 - value,
            )
            for donor, values in [("x", [0.1, 0.2, 0.8]), ("y", [0.4, 0.7, 0.9])]
            for comp in ["CD4", "CD8"]
            for i, value in enumerate(values)
        ]
    )


def test_sensitivity_summarizes_within_donor_before_across_donors():
    draws = synthetic_replicates()
    summary = replicate_summary(draws, ["x", "y"], 3)
    row = summary[summary.donor.eq("x") & summary.compartment.eq("CD4")].iloc[0]
    assert row.rfu_tv == 0.2
    assert row.rfu_tv_lo == pytest.approx(0.105)
    assert row.rfu_tv_hi == pytest.approx(0.77)
    across = summarize_donors(summary, ["rfu_tv"], ["compartment"])
    assert across.n_donors.tolist() == [2, 2]
    assert np.allclose(across["median"], 0.45)
    assert not np.isclose(draws[draws.compartment.eq("CD4")].d_group.median(), 0.45)
    with pytest.raises(ValueError, match="within donor"):
        summarize_donors(draws, ["d_group"], ["compartment"])


def test_missing_or_duplicated_technical_replicates_cannot_silently_plot():
    draws = synthetic_replicates()
    with pytest.raises(ValueError, match="coverage"):
        replicate_summary(draws.iloc[1:], ["x", "y"], 3)
    bad = draws.copy()
    bad.loc[1, "replicate"] = 0
    with pytest.raises(ValueError, match="Duplicate frozen"):
        replicate_summary(bad, ["x", "y"], 3)


def test_frozen_manifest_checks_unused_outputs_too(tmp_path):
    root = tmp_path / "frozen"
    root.mkdir()
    used, unused = root / "summary.json", root / "original.svg"
    used.write_text("{}")
    unused.write_text("original artifact")
    manifest = root / "completion.json"
    manifest.write_text(
        json.dumps(
            {
                "status": "complete",
                "outputs": {used.name: sha256(used), unused.name: sha256(unused)},
            }
        )
    )
    config = {"upstream": {"frozen": sha256(manifest)}}
    assert verify_upstream(tmp_path, config)["frozen"]["verified_outputs"] == 2
    unused.write_text("altered")
    with pytest.raises(ValueError, match="Upstream hash failed"):
        verify_upstream(tmp_path, config)
    manifest.write_text("{}")
    with pytest.raises(ValueError, match="completion hash changed"):
        verify_upstream(tmp_path, config)


def test_manifest_paths_cannot_escape_and_columns_must_be_documented(tmp_path):
    assert safe_child(tmp_path, "tables/A.tsv") == tmp_path / "tables/A.tsv"
    with pytest.raises(ValueError, match="escapes"):
        safe_child(tmp_path, "../outside")
    assert "persistent_rfus" in column_note("fraction_persistent_without_shared_clone")
    with pytest.raises(ValueError, match="Undocumented"):
        column_note("unexplained_new_endpoint")
