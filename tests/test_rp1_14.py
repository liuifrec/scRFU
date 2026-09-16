"""Entirely synthetic RP1-14 schema/alignment tests; no participant data."""

import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from manuscript.scripts.rp1_14_longitudinal import (
    analysis_pairs,
    empirical_reads,
    total_variation,
    validate_measurement_stage,
)
from manuscript.scripts.rp1_14_prepare import (
    assert_fixed_mapping,
    normalize_raw,
    receptor_features,
    reconcile_matrix,
    registry_from_metadata,
    restore_alignment,
)
from scrfu.tl import multiscale_repertoire_change


def source():
    # A non-C sequence in the middle causes the historical prefix/vector shift.
    return pd.DataFrame(
        {
            "aminoAcid": ["CAS", "AS", "CAT", "CAS"],
            "vMaxResolved": ["V1"] * 4,
            "jMaxResolved": ["J1"] * 4,
            "vGeneName": ["V1"] * 4,
            "jGeneName": ["J1"] * 4,
            "nucleotide": ["TGTGCTTCT", "GCTTCT", "TGTGCTACT", "TGTGCTTCT"],
            "vIndex": [0] * 4,
            "cdr3Length": [9, 6, 9, 9],
            "read_count": [40, 30, 20, 10],
            "sequenceStatus": ["In"] * 4,
        }
    )


def exported():
    return pd.DataFrame(
        {
            "cdr3_aa": ["CAS", "AS", "CAT"],
            "trbv": ["V1"] * 3,
            "rfu": [1, 2, 1],
            "max_cor": [0.8, 0.5, 0.8],
            "pass_thr": [True, False, True],
            "freq": [np.nan] * 3,
        }
    )


def test_registry_orders_by_age_and_preserves_missing_compartment_visit():
    metadata = pd.DataFrame(
        {
            "Sample ID": ["D901-CD4_1", "D901-CD4_3", "D901-CD8_3"],
            "CD4 or CD8": ["CD4", "CD4", "CD8"],
            "Age at collection": [50.0, 30.0, 30.0],
        }
    )
    result = registry_from_metadata(metadata)
    assert len(result) == 3
    assert result.loc[result.source_visit.eq("3"), "visit"].eq(1).all()
    assert result.loc[result.source_visit.eq("1"), "elapsed_years"].eq(20).all()
    assert set(result.donor) == {"RP01"}
    assert not ((result.compartment == "CD8") & (result.visit == 2)).any()


@pytest.mark.parametrize("change", ["label", "compartment", "age", "duplicate"])
def test_registry_rejects_mixed_or_conflicting_metadata(change):
    metadata = pd.DataFrame(
        {
            "Sample ID": ["D901-CD4_3", "D901-CD8_3"],
            "CD4 or CD8": ["CD4", "CD8"],
            "Age at collection": [30.0, 30.0],
        }
    )
    if change == "label":
        metadata.loc[0, "Sample ID"] = "RERF_144_unpublished_sample"
    elif change == "compartment":
        metadata.loc[0, "CD4 or CD8"] = "CD8"
    elif change == "age":
        metadata.loc[1, "Age at collection"] = 31
    else:
        metadata = pd.concat([metadata, metadata.iloc[[0]]])
    with pytest.raises(ValueError):
        registry_from_metadata(metadata)


@pytest.mark.parametrize("header", ["count", "count (reads)"])
def test_read_schema_and_denominator(header):
    raw = pd.DataFrame({header: [40, 30, 20, 10], "frequencyCount (%)": [40.0, 30.0, 20.0, 10.0]})
    result = normalize_raw(raw)
    assert result.read_count.sum() == 100
    repaired, _ = restore_alignment(source(), exported(), 4)
    assert repaired.read_count.sum() == 70  # no renormalization or imputed frequencies


def test_transformed_or_frequency_only_input_rejected():
    for raw in [
        pd.DataFrame({"frequencyCount (%)": [50, 50]}),
        pd.DataFrame({"count": [0.3, 0.7]}),
        pd.DataFrame({"count": [-1, 3]}),
        pd.DataFrame({"count": [5, 5], "frequencyCount (%)": [80, 20]}),
    ]:
        with pytest.raises(ValueError):
            normalize_raw(raw)


def test_alignment_repair_preserves_vectors_and_restores_fixed_map():
    saved = exported()
    repaired, audit = restore_alignment(source(), saved, 4)
    assert repaired.aminoAcid.tolist() == ["CAS", "CAT", "CAS"]
    assert repaired.rfu_numeric.tolist() == saved.rfu.tolist()
    assert repaired.max_cor.tolist() == saved.max_cor.tolist()
    assert repaired.source_row.tolist() == [0, 2, 3]
    assert audit["mislabelled_aa_rows_repaired"] == 2
    assert_fixed_mapping(receptor_features(repaired))


def test_unverified_repair_or_invalid_labels_rejected():
    for saved in [
        exported().assign(cdr3_aa=["CAT", "CAS", "CAS"]),
        exported().iloc[:2],
        exported().assign(rfu=[0, 2, 1]),
        exported().assign(pass_thr=True),
    ]:
        with pytest.raises(ValueError):
            restore_alignment(source(), saved, 4)


def test_nucleotide_cdr3_requires_translation_not_just_slice_length():
    features = receptor_features(source())
    assert features.clone.nunique() == 3
    with pytest.raises(ValueError, match="translation"):
        receptor_features(source().assign(aminoAcid="CAS"))


def test_fixed_mapping_conflicts_fail_before_distances():
    repaired, _ = restore_alignment(source(), exported(), 4)
    x = receptor_features(repaired)
    bad = x.copy()
    bad.loc[2, "rfu"] = "RFU2"
    with pytest.raises(ValueError, match="fixed assignment"):
        assert_fixed_mapping(bad)
    mapping = x.drop_duplicates("clone").set_index("clone").rfu
    a = x.groupby("clone").read_count.sum()
    b = a.iloc[::-1].copy()
    b.iloc[0] += 10
    result = multiscale_repertoire_change(a, b, mapping)
    assert result.summary["d_group"] <= result.summary["d_clone"] + 1e-12


def test_historical_matrix_is_assignment_row_mass_not_reads():
    aligned, _ = restore_alignment(source(), exported(), 4)
    values = np.zeros(5000)
    values[0], values[1] = 2 / 3 * 10000, 1 / 3 * 10000
    matrix = pd.DataFrame({"synthetic": values})
    assert reconcile_matrix(matrix, {"synthetic": aligned}, 10000) < 1e-10
    with pytest.raises(ValueError, match="transformed"):
        reconcile_matrix(np.log1p(matrix), {"synthetic": aligned}, 10000)
    with pytest.raises(ValueError, match="authorized"):
        reconcile_matrix(matrix.assign(unapproved=values), {"synthetic": aligned}, 10000)


def test_biological_inputs_ignored_without_hiding_schema():
    repo = Path(__file__).resolve().parents[1]
    paths = [
        "data/Data to be shared_RP P1-14/synthetic.xlsx",
        "data/RFU_out_RP_P1-14/synthetic.tsv",
        "data/schema.md",
    ]
    p = subprocess.run(
        ["git", "check-ignore", "--no-index", "--stdin"],
        cwd=repo,
        input="\n".join(paths) + "\n",
        text=True,
        capture_output=True,
        check=False,
    )
    assert p.returncode == 0
    assert set(p.stdout.splitlines()) == set(paths[:2])
    tracked = subprocess.check_output(["git", "ls-files", "data"], cwd=repo, text=True)
    assert not any(x in tracked for x in ["Data to be shared_RP P1-14", "RFU_out_RP_P1-14"])


def test_rp_count_units_and_controls_retain_frozen_development_settings():
    repo = Path(__file__).resolve().parents[1]
    rp = json.loads((repo / "manuscript/config/rp1_14_v1.json").read_text())
    dev = json.loads((repo / rp["frozen_development_config"]).read_text())
    assert rp["dataset_id"] == "RP1-14_published_authorized"
    assert rp["depth_cap_reads"] == dev["depth_cap_cells"]
    for key in ["depth_replicates", "random_group_replicates", "random_seed"]:
        assert rp[key] == dev[key]


def test_chronological_pairs_do_not_impute_missing_visit():
    coverage = pd.DataFrame(
        {
            "sample_id": ["synthetic_pre", "synthetic_late", "synthetic_other"],
            "donor": ["RP01", "RP01", "RP02"],
            "compartment": ["CD4"] * 3,
            "visit": [1, 3, 1],
            "elapsed_years": [0.0, 20.0, 0.0],
        }
    )
    result = analysis_pairs(coverage)
    assert len(result) == 1
    assert result.primary_pair.all() and result.interval_years.iloc[0] == 20
    assert result.visit_after.iloc[0] == 3


def test_read_resampling_uses_integer_counts_and_is_reproducible():
    counts = pd.Series([40, 10, 0], index=["a", "b", "c"])
    a = empirical_reads(counts, 20, np.random.default_rng(12))
    b = empirical_reads(counts, 20, np.random.default_rng(12))
    pd.testing.assert_series_equal(a, b)
    assert a.sum() == 20 and a.le(counts.loc[a.index]).all()
    pd.testing.assert_series_equal(
        empirical_reads(counts, 50, np.random.default_rng(12)), counts[counts > 0]
    )
    with pytest.raises(ValueError, match="integer"):
        empirical_reads(counts / 100, 1, np.random.default_rng(1))
    with pytest.raises(ValueError, match="exceeds"):
        empirical_reads(counts, 51, np.random.default_rng(1))


def test_full_clone_distance_retains_unmapped_mass_separately():
    a = pd.Series([8, 2], index=["mapped", "unmapped"])
    b = pd.Series([2, 8], index=["mapped", "unmapped"])
    assert total_variation(a, b) == pytest.approx(0.6)
    conditional = multiscale_repertoire_change(
        a, b, pd.Series({"mapped": "RFU1"}), unmapped="condition"
    )
    assert conditional.summary["d_clone"] == 0
    assert conditional.summary["coverage_before"] == 0.8
    assert conditional.summary["coverage_after"] == 0.2


def test_export_resumption_requires_a_verified_completed_measurement_stage(tmp_path):
    from scrfu.io import file_sha256

    identity = {"inputs": "synthetic", "settings": "frozen", "functions": "unchanged"}
    assert not validate_measurement_stage(tmp_path, identity)
    table = tmp_path / "synthetic.tsv"
    table.write_text("value\n1\n")
    marker = tmp_path / "measurement_stage.json"
    stage = {
        "status": "complete",
        "identity": identity,
        "outputs": {table.name: file_sha256(table)},
    }
    marker.write_text(json.dumps(stage))
    assert validate_measurement_stage(tmp_path, identity)
    with pytest.raises(ValueError, match="changed scientific"):
        validate_measurement_stage(tmp_path, {**identity, "settings": "different"})
    table.write_text("value\n2\n")
    with pytest.raises(ValueError, match="checksum"):
        validate_measurement_stage(tmp_path, identity)
