from __future__ import annotations

from types import SimpleNamespace

import pandas as pd
import pytest

from scrfu.adapters import (
    adapt_airr_dataframe,
    adapt_wells_tcr_ir,
    get_receptor_adapter,
    list_receptor_adapters,
)
from scrfu.pp import canonicalize_receptor_table, validate_receptor_table


def _minimal() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "input_row_id": ["row_00000000", "row_00000001"],
            "cell_id": ["c1", "c2"],
            "chain": ["TRB", "tra"],
            "cdr3aa": ["CASSA", "CAVR"],
            "v_call": ["TRBV1", "TRAV1"],
            "productive": ["true", 1],
            "source_adapter": ["synthetic", "synthetic"],
            "source_row_id": ["10", "11"],
        }
    )


def test_minimal_schema_normalization_and_source_order() -> None:
    normalized = canonicalize_receptor_table(_minimal())
    report = validate_receptor_table(normalized)

    assert normalized["input_row_id"].tolist() == ["row_00000000", "row_00000001"]
    assert normalized["chain"].tolist() == ["TRB", "TRA"]
    assert normalized["productive"].tolist() == [True, True]
    assert report["status"] == "ok"
    assert report["row_count"] == 2
    assert report["unique_cell_count"] == 2


@pytest.mark.parametrize(
    ("column", "value", "message"),
    [
        ("input_row_id", "row_00000000", "duplicates"),
        ("cell_id", None, "cell_id contains missing"),
        ("cdr3aa", None, "cdr3aa contains missing"),
    ],
)
def test_schema_actionable_errors(column: str, value: object, message: str) -> None:
    table = canonicalize_receptor_table(_minimal())
    table.loc[1, column] = value
    with pytest.raises(ValueError, match=message):
        validate_receptor_table(table)
    assert validate_receptor_table(table, strict=False)["status"] == "error"


def test_schema_optional_columns_are_not_required_and_source_duplicates_are_reported() -> None:
    table = canonicalize_receptor_table(_minimal())
    table["source_row_id"] = "same"
    report = validate_receptor_table(table)

    assert report["duplicated_source_row_ids"] == ["same"]
    assert "umi_count" not in table


def _airr() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "barcode": ["c1", "c1", "c2", "c2"],
            "locus": ["TRB", "TRB", "TRA", "TRB"],
            "junction_aa": ["CASSLOW", "CASSHIGH", "CAVR", "CASST"],
            "v_gene": ["TRBV1", "TRBV2", "TRAV1", "TRBV3"],
            "productive": [True, True, True, False],
            "umi_count": [1, 10, 5, 20],
        },
        index=["a", "b", "c", "d"],
    )


def test_airr_aliases_chain_filter_and_deterministic_primary_selection() -> None:
    result = adapt_airr_dataframe(_airr(), chain="TRB", primary_chain=True)

    assert result.receptors[["cell_id", "cdr3aa", "source_row_id"]].to_dict("records") == [
        {"cell_id": "c1", "cdr3aa": "CASSHIGH", "source_row_id": "b"}
    ]
    assert result.qc["rows_before_primary_selection"] == 2


def test_airr_all_chain_mode_without_counts_preserves_source_order() -> None:
    frame = _airr().drop(columns=["umi_count", "productive"])
    result = adapt_airr_dataframe(frame, chain=None, productive_only=False, primary_chain=False)

    assert result.receptors["source_row_id"].tolist() == ["a", "b", "c", "d"]
    assert result.receptors["chain"].tolist() == ["TRB", "TRB", "TRA", "TRB"]


def test_airr_metadata_remains_separate_and_barcode_mismatch_is_rejected() -> None:
    obs = pd.DataFrame({"donor": ["d1", "d2"]}, index=["c1", "c2"])
    result = adapt_airr_dataframe(_airr(), metadata=obs, metadata_columns=["donor"])
    assert result.cell_metadata.columns.tolist() == ["donor"]
    assert "donor" not in result.receptors

    with pytest.raises(ValueError, match="absent from cell metadata"):
        adapt_airr_dataframe(_airr(), metadata=obs.iloc[:1])


def test_wells_primary_priority_all_chain_mode_and_non_c_preservation() -> None:
    cells = pd.Index(["c1", "c2"])
    tcr = pd.DataFrame(
        {
            "TCR-IR_VDJ_1_locus": ["TRB", "TRB"],
            "TCR-IR_VDJ_1_junction_aa": ["CASS1", "ASS2"],
            "TCR-IR_VDJ_1_v_call": ["TRBV1", "TRBV2"],
            "TCR-IR_VDJ_1_productive": [True, True],
            "TCR-IR_VDJ_1_duplicate_count": [1, 5],
            "TCR-IR_VDJ_1_consensus_count": [9, 1],
            "TCR-IR_VDJ_2_locus": ["TRB", "TRA"],
            "TCR-IR_VDJ_2_junction_aa": ["CASS2", "CAV"],
            "TCR-IR_VDJ_2_v_call": ["TRBV3", "TRAV1"],
            "TCR-IR_VDJ_2_productive": [True, True],
            "TCR-IR_VDJ_2_duplicate_count": [2, 10],
            "TCR-IR_VDJ_2_consensus_count": [1, 1],
        },
        index=cells,
    )
    adata = SimpleNamespace(
        uns={"TCR_IR": tcr}, obsm={}, obs_names=cells, obs=pd.DataFrame(index=cells)
    )

    primary = adapt_wells_tcr_ir(adata)
    all_rows = adapt_wells_tcr_ir(adata, primary_chain=False)

    assert primary.receptors[["cell_id", "cdr3aa", "source_slot"]].to_dict("records") == [
        {"cell_id": "c1", "cdr3aa": "CASS2", "source_slot": "VDJ_2"},
        {"cell_id": "c2", "cdr3aa": "ASS2", "source_slot": "VDJ_1"},
    ]
    assert all_rows.receptors["cdr3aa"].tolist() == ["CASS1", "ASS2", "CASS2"]


def test_adapter_registry_is_explicit() -> None:
    assert list_receptor_adapters() == [
        "cellranger_vdj",
        "generic_airr_dataframe",
        "scirpy_airr",
        "wells_tcr_ir",
    ]
    assert get_receptor_adapter("wells") is adapt_wells_tcr_ir
    with pytest.raises(ValueError, match="Unknown receptor adapter"):
        get_receptor_adapter("guess")
