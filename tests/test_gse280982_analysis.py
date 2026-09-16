"""Synthetic external-registry/support tests; no biological downloads or R."""

import json

import numpy as np
import pandas as pd
import pytest

from manuscript.scripts.gse280982_analysis import (
    REPO,
    attach,
    check_cells,
    frozen_config,
    measurements,
    support_tables,
    verified_stage,
)
from manuscript.scripts.radiation_methods import json_save


def synthetic():
    rows, registry = [], []
    for tissue, visit, qualified in [
        ("tumor", 1, 100),
        ("tumor", 2, 99),
        ("tumor", 3, 101),
        ("blood", 1, 110),
    ]:
        sample = f"synthetic_{tissue}_{visit}"
        registry.append(
            dict(
                sample_id=sample,
                donor="synthetic",
                compartment=tissue,
                visit=visit,
                dataset_id="GSE280982",
                authorized=True,
                status="prepared_before_RFU_assignment",
                primary_TRB_cells_in_GEX=120,
                gex_barcodes=200,
            )
        )
        for i in range(120):
            rows.append(
                dict(
                    sample_id=sample,
                    donor="synthetic",
                    compartment=tissue,
                    visit=visit,
                    source_barcode=f"bc{i}",
                    cell_id=f"{sample}:bc{i}",
                    input_row_id=f"{sample}:row{i}",
                    dataset_id="GSE280982",
                    rfu_pass_threshold=i < qualified,
                    clone=f"clone{i % 6}",
                    rfu=f"RFU{i % 3 + 1}",
                    v_call="V",
                    j_call="J",
                    vj="V|J",
                    length_bin=2,
                )
            )
    return pd.DataFrame(rows), pd.DataFrame(registry)


def test_missing_visits_and_tissues_never_combined_or_zero_filled():
    frame, registry = synthetic()
    check_cells(frame, registry)
    visits, pairs = support_tables(frame, registry)
    assert len(pairs) == 6
    assert (
        visits.query("compartment == 'tumor' and visit == 2").status.item()
        == "below_frozen_support"
    )
    assert visits.query("compartment == 'blood' and visit == 2").qualified_cells.isna().all()
    good = pairs[pairs.status.eq("analyzed")]
    assert len(good) == 1
    assert (good.compartment.item(), good.visit_before.item(), good.visit_after.item()) == (
        "tumor",
        1,
        3,
    )
    with pytest.raises(ValueError, match="100-qualified"):
        support_tables(frame, registry, minimum=99)


def test_barcode_namespace_and_sample_metadata_enforced():
    frame, registry = synthetic()
    with pytest.raises(ValueError, match="namespace"):
        check_cells(frame.assign(cell_id=frame.cell_id.str.replace("synthetic", "other")), registry)
    with pytest.raises(ValueError, match="metadata"):
        check_cells(frame.assign(visit=1), registry)
    with pytest.raises(ValueError, match="unique"):
        check_cells(pd.concat([frame, frame.iloc[:1]]), registry)


def test_full_frozen_config_identity_not_selected_parameter_check(tmp_path):
    config = json.loads((REPO / "manuscript/config/radiation_methods_v1.json").read_text())
    path = tmp_path / "frozen_development_configuration.json"
    json_save(config, path)
    assert frozen_config(tmp_path)["reference"]["rfu_threshold"] == 0.6
    for key in ["depth_replicates", "min_cells_per_visit"]:
        changed = {**config, key: config[key] - 1}
        json_save(changed, path)
        with pytest.raises(ValueError, match="checksum"):
            frozen_config(tmp_path)


def test_nearest_cannot_rescue_unsupported_pair_and_cells_not_contig_reads(tmp_path):
    frame, registry = synthetic()
    _, pairs = support_tables(frame, registry)
    config = dict(
        random_seed=7,
        random_group_replicates=2,
        depth_replicates=2,
        depth_cap_cells=100,
        weights=["cell", "unique_clone"],
        rfu_detection_cell_counts=[1, 5],
    )
    measurements(frame.assign(read_count=10_000), pairs, config, tmp_path)
    result = pd.read_csv(tmp_path / "multiscale_pairs.tsv", sep="\t")
    assert len(result) == 12
    assert set(result.interval) == {"pre_to_6weeks"}
    assert set(result.compartment) == {"tumor"}
    assert set(result.policy) == {"threshold", "nearest"}
    assert (result.d_group <= result.d_clone + 1e-12).all()
    depth = pd.read_csv(tmp_path / "cell_depth_sensitivity.tsv", sep="\t")
    assert set(depth.cells_per_visit) == {100}
    before = (tmp_path / "fixed_group_controls.tsv").read_bytes()
    measurements(frame, pairs, config, tmp_path)
    assert before == (tmp_path / "fixed_group_controls.tsv").read_bytes()


def test_mapping_conflict_and_missing_assignment_rejected():
    frame = pd.DataFrame(
        dict(
            input_row_id=["a", "b"],
            junction=["nt", "nt"],
            v_call=["V", "V"],
            j_call=["J", "J"],
            cdr3aa=["CAS", "CAS"],
        )
    )
    rows = pd.DataFrame(
        dict(
            input_row_id=["a", "b"],
            rfu_label_nearest=["RFU1", "RFU2"],
            rfu_score=[0.7, 0.7],
            rfu_pass_threshold=[True, True],
            eligibility_status=["eligible", "eligible"],
        )
    )
    with pytest.raises(ValueError, match="varies"):
        attach(frame, rows)
    with pytest.raises(ValueError, match="coverage"):
        attach(frame, rows.iloc[:1])


def test_no_completion_from_evidence_counts_alone(tmp_path):
    json_save({"primary_TRB_cells": 100}, tmp_path / "evidence_counts.json")
    json_save({"status": "running", "fingerprint": "x"}, tmp_path / "completion.json")
    with pytest.raises(ValueError, match="not complete"):
        verified_stage(tmp_path)


def test_same_group_turnover_is_not_receptor_persistence(tmp_path):
    from manuscript.scripts.gse280982_analysis import persistence

    frame = pd.DataFrame(dict(visit=[1, 1, 3, 3], clone=["a", "a", "b", "b"], rfu=["RFU1"] * 4))
    row = persistence(frame, {"visit_before": 1, "visit_after": 3}, [1])[0]
    assert row["observed_before"] and row["observed_after"]
    assert row["observed_shared_receptors"] == 0
    assert np.isclose(row["dominant_fraction_after"], 1)
