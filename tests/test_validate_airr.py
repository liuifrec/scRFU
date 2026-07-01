from __future__ import annotations

import pandas as pd
import pytest

from scrfu.tl import validate_airr


class DummyAdata:
    def __init__(self, obs_names: list[str], airr: pd.DataFrame | None = None) -> None:
        self.obs_names = obs_names
        self.obs = pd.DataFrame(index=pd.Index(obs_names))
        self.obsm = {}
        if airr is not None:
            self.obsm["airr"] = airr


def test_validate_airr_valid_table() -> None:
    adata = DummyAdata(
        ["c1", "c2", "c3"],
        pd.DataFrame(
            {
                "barcode": ["c1", "c2", "c4"],
                "locus": ["TRB", "TRA", "TRB"],
                "cdr3": ["CASSA", "CAVR", "CASST"],
                "v": ["TRBV1", "TRAV1", "TRBV2"],
                "productive": ["true", "false", "yes"],
            }
        ),
    )

    result = validate_airr(adata)

    row = result.iloc[0]
    assert row["status"] == "ok"
    assert row["n_obs_cells"] == 3
    assert row["n_airr_rows"] == 3
    assert row["n_airr_cells"] == 3
    assert row["n_overlap_cells"] == 2
    assert row["barcode_overlap_rate"] == pytest.approx(2 / 3)
    assert row["has_cell_id_col"]
    assert row["has_chain_col"]
    assert row["has_cdr3_col"]
    assert row["has_v_col"]
    assert row["has_productive_col"]
    assert row["n_chain_rows"] == 2
    assert row["n_productive_chain_rows"] == 2


def test_validate_airr_missing_key_returns_status() -> None:
    result = validate_airr(DummyAdata(["c1"]))

    row = result.iloc[0]
    assert row["status"] == "missing_airr_key"
    assert row["n_obs_cells"] == 1
    assert row["n_airr_rows"] == 0


def test_validate_airr_missing_cdr3_column() -> None:
    adata = DummyAdata(
        ["c1"],
        pd.DataFrame({"cell_id": ["c1"], "chain": ["TRB"], "v_call": ["TRBV1"]}),
    )

    result = validate_airr(adata)

    row = result.iloc[0]
    assert row["status"] == "missing_required_columns"
    assert not row["has_cdr3_col"]
    assert row["has_v_col"]


def test_validate_airr_missing_v_column() -> None:
    adata = DummyAdata(
        ["c1"],
        pd.DataFrame({"cell_id": ["c1"], "chain": ["TRB"], "cdr3aa": ["CASSA"]}),
    )

    result = validate_airr(adata)

    row = result.iloc[0]
    assert row["status"] == "missing_required_columns"
    assert row["has_cdr3_col"]
    assert not row["has_v_col"]


def test_validate_airr_no_trb_rows() -> None:
    adata = DummyAdata(
        ["c1"],
        pd.DataFrame(
            {
                "cell_id": ["c1"],
                "chain": ["TRA"],
                "cdr3aa": ["CAVR"],
                "v_call": ["TRAV1"],
                "productive": [True],
            }
        ),
    )

    result = validate_airr(adata, chain="TRB")

    row = result.iloc[0]
    assert row["status"] == "no_chain_rows"
    assert row["n_chain_rows"] == 0


def test_validate_airr_barcode_mismatch() -> None:
    adata = DummyAdata(
        ["c1", "c2"],
        pd.DataFrame(
            {
                "cell_id": ["x1", "x2"],
                "chain": ["TRB", "TRB"],
                "cdr3aa": ["CASSA", "CASST"],
                "v_call": ["TRBV1", "TRBV2"],
            }
        ),
    )

    result = validate_airr(adata)

    row = result.iloc[0]
    assert row["status"] == "barcode_mismatch"
    assert row["n_overlap_cells"] == 0
    assert row["barcode_overlap_rate"] == 0.0


def test_validate_airr_strict_raises_for_missing_key() -> None:
    with pytest.raises(ValueError, match="missing obsm"):
        validate_airr(DummyAdata(["c1"]), strict=True)


def test_validate_airr_strict_raises_for_missing_required_columns() -> None:
    adata = DummyAdata(["c1"], pd.DataFrame({"cell_id": ["c1"]}))

    with pytest.raises(ValueError, match="missing required columns"):
        validate_airr(adata, strict=True)
