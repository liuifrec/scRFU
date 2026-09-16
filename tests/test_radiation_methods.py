from __future__ import annotations

import json

import pandas as pd
import pytest

from manuscript.scripts.radiation_methods import (
    assert_hash,
    assert_reference,
    completed,
    load_gse,
    reconcile_libraries,
    saved_rfu_labels,
    validate_boundary,
    wells_context,
)
from scrfu.io import file_sha256


def test_cohort_boundary_rejects_mixed_unapproved_and_unregistered_samples():
    frame = pd.DataFrame(
        {"dataset_id": ["RP1-14_published_authorized"], "sample_id": ["synthetic_allowed"]}
    )
    registry = frame.assign(authorized=True)
    validate_boundary(frame, registry)
    for dataset in ["RERF_144_unpublished", "mixed_RP_RERF", "unknown"]:
        with pytest.raises(ValueError, match="Unapproved or mixed"):
            validate_boundary(frame.assign(dataset_id=dataset), registry)
    with pytest.raises(ValueError, match="authorized"):
        validate_boundary(frame.assign(sample_id="not_in_subset"), registry)
    with pytest.raises(ValueError, match="authorization"):
        validate_boundary(frame.drop(columns="dataset_id"), registry)
    with pytest.raises(ValueError, match="authorized"):
        validate_boundary(frame, registry.assign(authorized=False))


def test_cache_requires_completion_hashes_and_settings(tmp_path):
    (tmp_path / "table.tsv").write_text("x\n1\n")
    assert not completed(tmp_path, "scientific_fingerprint")
    manifest = {
        "status": "running",
        "fingerprint": "scientific_fingerprint",
        "outputs": {"table.tsv": file_sha256(tmp_path / "table.tsv")},
    }
    status = tmp_path / "completion.json"
    status.write_text(json.dumps(manifest))
    assert not completed(tmp_path, "scientific_fingerprint")
    manifest["status"] = "complete"
    status.write_text(json.dumps(manifest))
    assert completed(tmp_path, "scientific_fingerprint")
    with pytest.raises(ValueError, match="different scientific"):
        completed(tmp_path, "changed_settings")
    (tmp_path / "table.tsv").write_text("x\n2\n")
    with pytest.raises(ValueError, match="changed"):
        completed(tmp_path, "scientific_fingerprint")


def test_checksum_and_reference_backend_must_match(tmp_path):
    path = tmp_path / "synthetic.tsv"
    path.write_text("x\n1\n")
    assert_hash(path, file_sha256(path))
    with pytest.raises(ValueError, match="checksum"):
        assert_hash(path, "0" * 64)
    manifest = {
        "backend_mode": "standard",
        "rfu_threshold": 0.6,
        "km5000_rdata_sha256": "synthetic",
        "chunk_count": 2,
        "completed_chunk_count": 2,
        "failed_chunk_count": 0,
    }
    expected = {k: manifest[k] for k in ["backend_mode", "rfu_threshold", "km5000_rdata_sha256"]}
    assert_reference(manifest, expected)
    with pytest.raises(ValueError, match="incompatibility"):
        assert_reference({**manifest, "backend_mode": "map_aware"}, expected)
    with pytest.raises(ValueError, match="incomplete"):
        assert_reference({**manifest, "completed_chunk_count": 1}, expected)


def test_saved_labels_are_not_zero_based_ids_or_renumbered():
    labels = pd.Series(["RFU1", "RFU5000"], dtype="string")
    pd.testing.assert_series_equal(saved_rfu_labels(labels), labels)
    with pytest.raises(ValueError, match="explicit saved"):
        saved_rfu_labels(pd.Series([0, 4999]))


def test_gse_metadata_conflict_rejected_before_analysis(tmp_path):
    # Only source layout is reproduced: these rows are wholly synthetic.
    raw = pd.DataFrame(
        {"Unnamed: 0": ["b1_x"], "patient": ["synthetic"], "state": ["before"], "method": ["SBRT"]}
    )
    meta = pd.DataFrame(
        {"Unnamed: 0": ["b1_x"], "patient": ["P001"], "state": ["Post"], "method": ["SBRT"]}
    )
    raw.to_csv(tmp_path / "GSE190905_TCR_data.csv.gz", index=False)
    meta.to_csv(tmp_path / "GSE190905_meta_data.csv.gz", index=False)
    (tmp_path / "rfu").mkdir()
    pd.DataFrame({"cell_id": ["b1_x"]}).to_csv(
        tmp_path / "rfu/rfu_results_per_row.tsv.gz", sep="\t", index=False
    )
    with pytest.raises(ValueError, match="treatment/timepoint"):
        load_gse(tmp_path, tmp_path / "not_needed.gz", {}, tmp_path)


def test_pooled_library_names_cannot_be_assumed_to_match():
    observations = pd.DataFrame(
        {
            "library": ["b4", "b4", "b5"],
            "patient_rna": ["synthetic_A", "synthetic_B", "synthetic_A"],
            "state_rna": ["Post", "Post", "Pre"],
        }
    )
    libraries = pd.concat(
        [
            observations.assign(
                library=["b8", "b8", "b1"],
                assay=assay,
                accession=["sample8_" + assay, "sample8_" + assay, "sample1_" + assay],
            )
            for assay in ("RNA", "TCR")
        ],
        ignore_index=True,
    )
    result = reconcile_libraries(observations, libraries)
    assert set(result.loc[result.processed_library.eq("b4"), "geo_library"]) == {"b8"}
    assert set(result.loc[result.processed_library.eq("b5"), "geo_library"]) == {"b1"}
    with pytest.raises(ValueError, match="unsupported"):
        reconcile_libraries(observations.assign(patient_rna="inconsistent"), libraries)
    with pytest.raises(ValueError, match="bijective"):
        reconcile_libraries(
            pd.concat([observations, observations.iloc[:2].assign(library="extra")]), libraries
        )


def test_wells_coverage_retains_atlas_donor_without_eligible_receptor(tmp_path):
    data = tmp_path / "full_run"
    data.mkdir()
    meta = pd.DataFrame(
        {
            "cell_id": ["s1", "s2", "s3"],
            "donor_id": ["synthetic_A", "synthetic_B", "synthetic_C"],
            "tissue": ["blood"] * 3,
            "cell_type": ["synthetic_state"] * 3,
            "library_id": ["lib1", "lib2", "lib3"],
        }
    )
    meta.to_csv(data / "obs_metadata.tsv.gz", sep="\t", index=False)
    rfu = pd.DataFrame(
        {
            "cell_id": ["s1", "s2"],
            "cdr3aa": ["CASS", "CAST"],
            "v_call": ["TRBV1", "TRBV2"],
            "rfu_pass_threshold": [True, True],
            "rfu_label_nearest": ["RFU1", "RFU2"],
        }
    )
    rfu.to_csv(data / "rfu_results_per_row.tsv.gz", sep="\t", index=False)
    summary = wells_context(tmp_path, {"wells_min_cells": 1, "wells_min_donors": 1}, tmp_path)
    coverage = pd.read_csv(tmp_path / "wells_donor_tissue_coverage.tsv", sep="\t")
    assert summary["atlas_donors"] == 3 and summary["donors"] == 2
    assert coverage.atlas_cells.sum() == 3
    assert len(coverage) == 3 and coverage.primary_trb_cells.eq(0).sum() == 1
