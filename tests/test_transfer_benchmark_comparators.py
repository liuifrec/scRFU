from __future__ import annotations

import pandas as pd
import pytest

import scrfu


def _hash(character: str) -> str:
    return character * 64


def _reference() -> scrfu.tl.FrozenRFUReference:
    return scrfu.tl.FrozenRFUReference.create(
        rfu_r_sha256=_hash("a"),
        embedding_sha256=_hash("b"),
        centroid_atlas_sha256=_hash("c"),
        threshold=0.6,
        eligibility_rule="starts_with_C",
        assignment_mode="standard",
        receptor_chain="TRB",
        receptor_model="original_rfu_trb",
        reference_label="test-reference",
        creation_timestamp="2026-01-01T00:00:00+00:00",
    )


def _rows() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "input_row_id": ["r1", "r2", "r3", "r4"],
            "cell_id": ["c1", "c2", "c3", "c4"],
            "sample": ["s1", "s1", "s2", "s2"],
            "cdr3aa": ["CASSF", "CASRF", "CASSF", "CATRF"],
            "v_call": ["TRBV1", "TRBV2", "TRBV1", "TRBV3"],
            "j_call": ["TRBJ1", "TRBJ1", "TRBJ1", "TRBJ2"],
            "clonotype_id": ["cl1", "cl2", "cl1", "cl3"],
            "rfu_label": ["RFU1", "RFU2", "RFU1", "RFU3"],
            "rfu_score": [0.8, 0.5, 0.9, 0.7],
            "reference_coverage_status": [
                "covered",
                "below_threshold",
                "covered",
                "covered",
            ],
            "pass_thr": [True, False, True, True],
        }
    )


def test_frozen_reference_identifier_validation() -> None:
    reference = _reference()
    assert reference.immutable_reference_id.startswith("scrfu-ref-")
    assert scrfu.tl.validate_frozen_reference(reference) is reference
    changed = reference.to_dict()
    changed["threshold"] = 0.7
    with pytest.raises(ValueError, match="identifier"):
        scrfu.tl.validate_frozen_reference(changed)
    changed["immutable_reference_id"] = reference.immutable_reference_id
    changed["rfu_r_sha256"] = "not-a-hash"
    with pytest.raises(ValueError, match="SHA256"):
        scrfu.tl.validate_frozen_reference(changed)


def test_transfer_summary_does_not_refit_and_checks_reference() -> None:
    reference = _reference()
    result = scrfu.tl.transfer_cohort(
        _rows(),
        reference,
        cohort_label="synthetic",
        sample_key="sample",
        observed_reference_id=reference.immutable_reference_id,
    )
    assert result.sample_matrix.shape == (2, 3)
    assert result.provenance["no_reference_refitting"] is True
    assert result.provenance["target_reference_verified"] is True
    assert result.coverage["input_rows"].sum() == 4
    assert result.rfu_summary["unique_sequence_count"].sum() == 3
    pd.testing.assert_frame_equal(result.source_table, _rows())
    with pytest.raises(ValueError, match="not generated"):
        scrfu.tl.transfer_cohort(
            _rows(),
            reference,
            cohort_label="synthetic",
            sample_key="sample",
            observed_reference_id="wrong",
        )
    incompatible = _rows().assign(chain="IGH")
    with pytest.raises(ValueError, match="incompatible"):
        scrfu.tl.transfer_cohort(
            incompatible,
            reference,
            cohort_label="synthetic",
            sample_key="sample",
        )


def test_explicit_metadata_harmonization_and_heldout_manifest() -> None:
    metadata = pd.DataFrame(
        {"subject": ["d1", "d2"], "sex_raw": ["female", "unknown"], "unused": [1, 2]}
    )
    result = scrfu.tl.harmonize_cohort_metadata(
        metadata,
        cohort_label="cohort-a",
        field_mapping={"donor": "subject", "sex": "sex_raw"},
        value_mappings={"sex": {"female": "F"}},
        strict=False,
    )
    assert result.metadata["sex"].tolist()[0] == "F"
    assert pd.isna(result.metadata["sex"].tolist()[1])
    assert result.unmapped_values["source_value"].tolist() == ["unknown"]
    assert result.rules["no_inferred_mappings"] is True
    with pytest.raises(ValueError, match="Unmapped values"):
        scrfu.tl.harmonize_cohort_metadata(
            metadata,
            cohort_label="cohort-a",
            field_mapping={"sex": "sex_raw"},
            value_mappings={"sex": {"female": "F"}},
        )

    manifest = scrfu.tl.create_heldout_validation_manifest(
        development_cohorts=["dev"],
        validation_cohorts=["validation"],
        held_out_cohort="heldout",
        reference=_reference(),
        frozen_parameters={"threshold": 0.6},
        data_hashes={
            "dev": _hash("d"),
            "validation": _hash("f"),
            "heldout": _hash("e"),
        },
        evaluation_metrics=["coverage", "donor_retrieval"],
    )
    assert manifest.held_out_cohort == "heldout"
    with pytest.raises(ValueError, match="overlap"):
        scrfu.tl.create_heldout_validation_manifest(
            development_cohorts=["same"],
            validation_cohorts=[],
            held_out_cohort="same",
            reference=_reference(),
            frozen_parameters={},
            data_hashes={},
            evaluation_metrics=["coverage"],
        )


def test_deterministic_subsampling_multinomial_and_stability() -> None:
    rows = _rows()
    first = scrfu.tl.deterministic_subsample(rows, unit="sequence", n=2, random_state=3)
    second = scrfu.tl.deterministic_subsample(rows, unit="sequence", n=2, random_state=3)
    pd.testing.assert_frame_equal(first, second)
    assert first["cdr3aa"].nunique() == 2

    matrix = pd.DataFrame([[3, 1], [0, 4]], index=["s1", "s2"], columns=["a", "b"])
    resampled_a = scrfu.tl.multinomial_abundance_resample(matrix, depth=20, random_state=7)
    resampled_b = scrfu.tl.multinomial_abundance_resample(matrix, depth=20, random_state=7)
    pd.testing.assert_frame_equal(resampled_a, resampled_b)
    assert resampled_a.sum(axis=1).tolist() == [20, 20]
    assert resampled_a.attrs["resampling"]["physical_read_downsampling"] is False

    stability = scrfu.tl.benchmark_representation_stability(matrix, matrix, top_k=1)
    values = stability.metrics.set_index(["sample", "metric"])["value"]
    assert values.loc[("s1", "cosine")] == pytest.approx(1)
    assert values.loc[("s1", "mean_absolute_error")] == 0
    assert values.loc[("s2", "top_k_overlap")] == 1

    sensitivity = scrfu.tl.threshold_sensitivity(_rows(), [0.5, 0.8])
    assert sensitivity["threshold_pass_count"].tolist() == [4, 2]
    shuffled_a = scrfu.tl.shuffle_input_order(rows, random_state=9)
    shuffled_b = scrfu.tl.shuffle_input_order(rows, random_state=9)
    pd.testing.assert_frame_equal(shuffled_a, shuffled_b)
    assert shuffled_a.index.tolist() != rows.index.tolist()
    folds = scrfu.tl.donor_leave_one_out(
        pd.DataFrame({"donor": ["d1", "d1", "d2"], "visit": [0, 1, 0]}),
        donor_key="donor",
    )
    assert [donor for donor, _ in folds] == ["d1", "d2"]
    assert folds[0][1]["donor"].tolist() == ["d2"]


@pytest.mark.parametrize(
    "method",
    [
        "exact_cdr3",
        "clonotype",
        "v_gene",
        "j_gene",
        "cdr3_length",
        "edit_distance",
        "shannon",
        "simpson",
        "dominant_clone_fraction",
    ],
)
def test_comparator_representations_share_sample_interface(method: str) -> None:
    result = scrfu.tl.repertoire_representation(_rows(), sample_key="sample", method=method)
    assert result.matrix.index.tolist() == ["s1", "s2"]
    assert result.parameters["unit"] == "biological_sample"
    assert result.matrix.shape[1] >= 1
