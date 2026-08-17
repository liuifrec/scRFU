from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import scrfu


def test_reference_coverage_hand_calculated_and_grouped() -> None:
    frame = pd.DataFrame(
        {
            "sample": ["s1", "s1", "s1", "s2", "s2"],
            "cell_id": ["c1", "c2", "c3", "c4", "c5"],
            "unique_sequence_id": ["q1", "q2", "q3", "q4", pd.NA],
            "cdr3aa": ["CASSF", "CASRF", "BAD", "CATRF", pd.NA],
            "reference_coverage_status": [
                "covered",
                "below_threshold",
                "ineligible_sequence",
                "upstream_unassigned",
                "missing_sequence",
            ],
            "rfu_score": [0.8, 0.4, np.nan, np.nan, np.nan],
        }
    )
    result = scrfu.tl.reference_coverage(frame, groupby="sample")
    first = result.set_index("sample").loc["s1"]
    assert first["input_rows"] == 3
    assert first["eligible_rows"] == 2
    assert first["assigned_rows"] == 2
    assert first["threshold_pass_rows"] == 1
    assert first["threshold_pass_fraction"] == pytest.approx(0.5)
    assert first["unique_sequence_coverage"] == pytest.approx(0.5)
    assert first["cell_weighted_coverage"] == pytest.approx(0.5)
    second = result.set_index("sample").loc["s2"]
    assert second["upstream_unassigned_rows"] == 1
    assert second["missing_sequence_rows"] == 1


def test_reference_coverage_legacy_fields_and_invalid_status() -> None:
    legacy = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3"],
            "cdr3aa": ["CASSF", "CASRF", pd.NA],
            "eligibility_status": ["eligible", "eligible", "missing_cdr3"],
            "rfu_id": [1, 2, pd.NA],
            "pass_thr": [True, False, False],
            "rfu_score": [0.8, 0.4, np.nan],
        }
    )
    row = scrfu.tl.reference_coverage(legacy).iloc[0]
    assert row["threshold_pass_rows"] == 1
    assert row["below_threshold_rows"] == 1
    assert row["missing_sequence_rows"] == 1
    legacy["reference_coverage_status"] = ["covered", "invented", "missing_sequence"]
    with pytest.raises(ValueError, match="Unknown reference coverage"):
        scrfu.tl.reference_coverage(legacy)


def test_mudata_like_routing_preserves_modalities_and_checks_alignment() -> None:
    airr = pd.DataFrame(
        {
            "cell_id": ["c1", "c2"],
            "chain": ["TRB", "TRB"],
            "junction_aa": ["CASSF", "CASRF"],
            "v_call": ["TRBV1", "TRBV2"],
            "productive": [True, True],
        }
    )
    receptor_modality = SimpleNamespace(
        obsm={"airr": airr}, obs=pd.DataFrame(index=pd.Index(["c1", "c2"]))
    )
    metadata = pd.DataFrame({"donor": ["d1", "d2"]}, index=["c1", "c2"])
    metadata_modality = SimpleNamespace(obs=metadata)
    mdata = SimpleNamespace(mod={"airr": receptor_modality, "rna": metadata_modality})
    original_keys = tuple(mdata.mod)

    result = scrfu.adapters.prepare_receptors(
        mdata,
        adapter="scirpy_airr",
        modality="airr",
        metadata_modality="rna",
        metadata_columns=["donor"],
    )
    assert result.receptors["cell_id"].tolist() == ["c1", "c2"]
    assert result.cell_metadata["donor"].tolist() == ["d1", "d2"]
    assert result.provenance["mudata_modality"] == "airr"
    assert tuple(mdata.mod) == original_keys

    with pytest.raises(KeyError, match="modality 'protein'"):
        scrfu.adapters.prepare_receptors(mdata, adapter="scirpy_airr", modality="protein")
    unaligned = SimpleNamespace(
        mod={
            "airr": receptor_modality,
            "rna": SimpleNamespace(obs=pd.DataFrame(index=["c1"])),
        }
    )
    with pytest.raises(ValueError, match="absent"):
        scrfu.adapters.prepare_receptors(
            unaligned,
            adapter="scirpy_airr",
            modality="airr",
            metadata_modality="rna",
        )
