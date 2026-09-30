"""Scientific invariants for data-only figure tables, using synthetic inputs."""

import json

import numpy as np
import pandas as pd
import pytest

from manuscript.scripts import cross_application_reusable_figures as pack


def test_external_missing_intervals_are_not_zero_and_unexpected_estimates_fail():
    support = pd.DataFrame(
        {
            "donor": ["a", "b"],
            "compartment": ["tumor"] * 2,
            "interval": ["P-R"] * 2,
            "status": ["analyzed", "unavailable_required_visit"],
        }
    )
    distance = (
        support.iloc[:1]
        .drop(columns="status")
        .assign(d_clone=0.8, d_group=0.6, aggregation_cancellation=0.2)
    )
    result = pack.external_pairs(support, distance)
    assert result.loc[1, pack.DIST].isna().all()
    assert result.loc[0, "d_group"] == 0.6
    invented = pd.concat([distance, distance.assign(donor="b")], ignore_index=True)
    with pytest.raises(ValueError, match="Unavailable interval"):
        pack.external_pairs(support, invented)
    with pytest.raises(ValueError, match="Duplicate external pair"):
        pack.external_pairs(support, pd.concat([distance, distance]))


def test_wells_absent_combinations_and_observed_zero_remain_distinct():
    source = pd.DataFrame(
        {
            "donor": ["a", "a", "b"],
            "tissue": ["blood", "lung", "blood"],
            "atlas_cells": [20, 5, 15],
            "primary_trb_cells": [10, 0, 10],
            "threshold_cells": [8, 0, 5],
            "threshold_fraction_of_primary_trb": [0.8, np.nan, 0.5],
        }
    )
    result = pack.wells_grid(source).set_index(["donor", "tissue"])
    assert result.loc[("a", "lung"), "coverage_state"] == "atlas_zero_primary_TRB"
    assert result.loc[("a", "lung"), "atlas_cells"] == 5
    assert result.loc[("b", "lung"), "coverage_state"] == "no_atlas_sample"
    assert np.isnan(result.loc[("b", "lung"), "atlas_cells"])
    source.loc[1, "threshold_fraction_of_primary_trb"] = 0
    with pytest.raises(ValueError, match="undefined"):
        pack.wells_grid(source)


def sensitivity_inputs():
    donors = [f"d{i}" for i in range(6)]
    primary = pd.DataFrame({"donor": donors, "d_group": np.arange(6) / 10})
    unique = primary.assign(
        policy="threshold",
        weighting="unique_clone",
        grouping="RFU",
        subset="all_T",
        status="valid",
        d_clone=lambda x: x.d_group + 0.2,
        aggregation_cancellation=0.2,
    )
    depth = pd.DataFrame(
        [
            {
                "donor": donor,
                "replicate": replicate,
                "cells_per_visit": 404 if i == 0 else 500,
                "d_group": (replicate / 100 + i / 20),
            }
            for i, donor in enumerate(donors)
            for replicate in range(50)
        ]
    )
    dominance = primary.assign(retained_cell_fraction_before=0.6, retained_cell_fraction_after=0.4)
    return primary, unique, depth, dominance


def test_sampling_summaries_preserve_donor_and_technical_replicate_units():
    tables = sensitivity_inputs()
    result = pack.sensitivity_table(*tables)
    draws = result.query("condition == 'depth'").set_index("donor")
    assert len(draws) == 6
    assert draws.loc["d0", "rfu_tv"] == pytest.approx(0.245)
    assert draws.loc["d0", "low"] == pytest.approx(0.01225)
    assert draws.loc["d0", "high"] == pytest.approx(0.47775)
    assert draws.loc["d0", "cells_per_visit"] == 404
    assert (draws.replicates == 50).all()
    with pytest.raises(ValueError, match="replicate coverage"):
        pack.sensitivity_table(tables[0], tables[1], tables[2].iloc[1:], tables[3])


def test_changed_weighting_cancellation_is_rejected():
    _, frame, _, _ = sensitivity_inputs()
    frame.loc[0, "aggregation_cancellation"] = 0.15
    with pytest.raises(ValueError, match="Cancellation identity"):
        pack.primary(frame, "unique_clone")


def regulatory_inputs():
    tier = "A_same_variant_conditional_layers"
    candidates = pd.DataFrame(
        [
            {"variant": f"v{i}", "rfu_label": f"r{j}", "evidence_tier": tier}
            for i, count in enumerate([36, 5, 5, 5, 5, 6])
            for j in range(count)
        ]
    )
    targets = pd.DataFrame(
        [
            {
                "variant": f"v{i}",
                "rfu_label": f"r{rfu}",
                "layer": layer,
                "target": f"{layer}_{i}_{target}",
                "conditional_rank": 1,
                "target_class": ("non_TCR_gene" if i == 5 else "direct_TCR_gene")
                if layer == "eqtl"
                else np.nan,
                "evidence_tier": tier,
                "conditional_beta": 0.1,
            }
            for i, peaks in enumerate([2, 1, 2, 3, 1, 1])
            for rfu in range(2)
            for layer, count in [("eqtl", 1), ("caqtl", peaks)]
            for target in range(count)
        ]
    )
    return candidates, targets


def test_regulatory_records_deduplicate_rfu_joins_and_preserve_target_classes():
    candidates, targets = regulatory_inputs()
    matrix, records = pack.regulatory_matrix(candidates, targets)
    assert len(matrix) == 6 and len(records) == 16
    assert matrix.caqtl_records.sum() == 10
    assert matrix.target_class.tolist() == ["direct_TCR_gene"] * 5 + ["non_TCR_gene"]
    conflicting = pd.concat([targets, targets.iloc[:1].assign(conditional_beta=99)])
    with pytest.raises(ValueError, match="Conflicting molecular records"):
        pack.regulatory_matrix(candidates, conflicting)


def prediction_inputs():
    data = {
        "aliases": pd.DataFrame(
            {
                "patient_tcr": [f"source{i}" for i in range(6)],
                "donor": [f"alias{i}" for i in range(6)],
            }
        )
    }
    for design, mean, delta in [("ordinary", 1.850996, 0.004216), ("purged", 1.852679, 0.004317)]:
        rows, metrics = [], []
        for i in range(6):
            for weight in ["clone", "cell"]:
                baseline = mean + (i - 2.5) / 100
                extended = baseline + delta
                rows.append(
                    {
                        "held_out_donor": f"source{i}",
                        "weighting": weight,
                        "baseline_log_loss": baseline,
                        "extended_log_loss": extended,
                        "delta_log_loss": delta,
                    }
                )
                for model, value in [
                    ("TRBV_TRBJ_length", baseline),
                    ("TRBV_TRBJ_length_RFU", extended),
                ]:
                    metrics.append(
                        {
                            "held_out_donor": f"source{i}",
                            "weighting": weight,
                            "model": model,
                            "log_loss": value,
                        }
                    )
        data[f"{design}_deltas"] = pd.DataFrame(rows).sample(frac=1, random_state=9)
        data[f"{design}_metrics"] = pd.DataFrame(metrics).sample(frac=1, random_state=10)
    return data


def test_prediction_uses_saved_fold_deltas_and_existing_tcr_aliases():
    data = prediction_inputs()
    rows, summary = pack.prediction_tables(data)
    assert set(rows.donor) == {f"alias{i}" for i in range(6)}
    assert len(rows) == 24 and len(summary) == 2
    assert np.allclose(summary.delta_log_loss, [0.004216, 0.004317])
    assert summary.donors.eq(6).all()
    data["ordinary_deltas"]["delta_log_loss"] *= -1
    with pytest.raises(ValueError, match="delta sign"):
        pack.prediction_tables(data)


def test_missing_prediction_alias_fails_closed():
    data = prediction_inputs()
    data["aliases"] = data["aliases"].iloc[:5]
    with pytest.raises(ValueError, match="donor alias"):
        pack.prediction_tables(data)


def test_all_upstream_outputs_are_verified_even_when_not_plotted(tmp_path, monkeypatch):
    root = tmp_path / "completed"
    root.mkdir()
    source = root / "table.tsv"
    source.write_text("value\n1\n")
    protected = root / "protected.svg"
    protected.write_text("original protected figure")
    completion = root / "completion.json"
    completion.write_text(
        json.dumps(
            {"status": "complete", "outputs": {p.name: pack.sha256(p) for p in [source, protected]}}
        )
    )
    config = {
        "upstream": {
            "fixture": {
                "workspace": "methods",
                "directory": "completed",
                "sha256": pack.sha256(completion),
            }
        }
    }
    monkeypatch.setattr(pack, "SOURCES", {"table": ("fixture", "table.tsv")})
    data, _, verified = pack.verify_sources({"methods": tmp_path}, config)
    assert data["table"].value.tolist() == [1]
    assert verified["fixture"]["verified_outputs"] == 2
    protected.write_text("changed figure")
    with pytest.raises(ValueError, match="Upstream hash failed"):
        pack.verify_sources({"methods": tmp_path}, config)
