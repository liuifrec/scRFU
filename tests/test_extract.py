import pandas as pd
import pytest


class DummyAdata:
    def __init__(
        self,
        airr_df: pd.DataFrame | None = None,
        obs_names=(),
        *,
        obsm: dict | None = None,
        uns: dict | None = None,
    ):
        self.obsm = ({"airr": airr_df} if airr_df is not None else {}) | (obsm or {})
        self.uns = uns or {}
        self.obs_names = obs_names


def test_extract_trb_features_basic():
    from scrfu.extract import extract_trb_features

    airr = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3"],
            "chain": ["TRB", "TRA", "TRB"],
            "cdr3aa": ["CASSA", "CAVR", "CASST"],
            "v_call": ["TRBV1", "TRAV1", "TRBV2"],
            "productive": [True, True, True],
        }
    )
    adata = DummyAdata(airr, obs_names=["c1", "c2", "c3"])

    out = extract_trb_features(adata, airr_key="airr", chain="TRB")
    assert set(out.columns) == {"cell_id", "cdr3aa", "trbv"}
    assert list(out["cell_id"]) == ["c1", "c3"]
    assert list(out["trbv"]) == ["TRBV1", "TRBV2"]


def test_extract_missing_key():
    from scrfu.extract import extract_trb_features

    class A:
        obsm = {}

    with pytest.raises(KeyError):
        extract_trb_features(A(), airr_key="airr", chain="TRB")


def test_extract_wells_tcr_ir_from_uns_preserves_verified_cell_alignment():
    from scrfu.extract import extract_trb_features

    wells = pd.DataFrame(
        {
            "TCR-IR_VDJ_1_locus": ["TRB", "TRA", "TRB", "TRB", "TRB", "TRB"],
            "TCR-IR_VDJ_1_junction_aa": [
                " CASSA ",
                "CAVR",
                None,
                "ASSQ",
                "CASST",
                "CASSN",
            ],
            "TCR-IR_VDJ_1_v_call": ["TRBV1", "TRAV1", "TRBV2", "TRBV3", "", "nan"],
            "TCR-IR_VDJ_1_productive": ["True", "True", "True", "True", "True", "True"],
        },
        index=["c1", "c2", "c3", "c4", "c5", "c6"],
    )
    original = wells.copy(deep=True)
    adata = DummyAdata(obs_names=wells.index, uns={"TCR_IR": wells})

    out = extract_trb_features(adata, airr_key="TCR_IR")

    expected = pd.DataFrame({"cell_id": ["c1"], "cdr3aa": ["CASSA"], "trbv": ["TRBV1"]})
    pd.testing.assert_frame_equal(out, expected)
    pd.testing.assert_frame_equal(adata.uns["TCR_IR"], original)


def test_extract_wells_tcr_ir_mapping_uses_obs_order_for_default_index():
    from scrfu.extract import extract_wells_tcr_ir_features

    wells = {
        "IR_VDJ_1_locus": ["TRB", "TRB"],
        "IR_VDJ_1_junction_aa": ["CASSA", "CASST"],
        "IR_VDJ_1_v_call": ["TRBV1", "TRBV2"],
    }
    adata = DummyAdata(obs_names=["cell-b", "cell-a"], uns={"TCR_IR": wells})

    out = extract_wells_tcr_ir_features(adata)

    assert out.to_dict("records") == [
        {"cell_id": "cell-b", "cdr3aa": "CASSA", "trbv": "TRBV1"},
        {"cell_id": "cell-a", "cdr3aa": "CASST", "trbv": "TRBV2"},
    ]


def test_extract_wells_tcr_ir_supports_public_obsm_representation():
    from scrfu.extract import extract_trb_features

    wells = pd.DataFrame(
        {
            "TCR-IR_VDJ_1_junction_aa": ["CASSA"],
            "TCR-IR_VDJ_1_v_call": ["TRBV1"],
        },
        index=["c1"],
    )
    adata = DummyAdata(obs_names=["c1"], obsm={"TCR_IR": wells})

    out = extract_trb_features(adata, airr_key="TCR_IR")

    assert out.to_dict("records") == [{"cell_id": "c1", "cdr3aa": "CASSA", "trbv": "TRBV1"}]


def test_extract_wells_tcr_ir_prefers_uns_over_obsm():
    from scrfu.extract import extract_wells_tcr_ir_features

    def table(cdr3: str) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "TCR-IR_VDJ_1_junction_aa": [cdr3],
                "TCR-IR_VDJ_1_v_call": ["TRBV1"],
            },
            index=["c1"],
        )

    adata = DummyAdata(
        obs_names=["c1"],
        uns={"TCR_IR": table("CASS_UNS")},
        obsm={"TCR_IR": table("CASS_OBSM")},
    )

    out = extract_wells_tcr_ir_features(adata)

    assert out.loc[0, "cdr3aa"] == "CASS_UNS"


def test_extract_wells_tcr_ir_rejects_ambiguous_row_alignment():
    from scrfu.extract import extract_wells_tcr_ir_features

    wells = pd.DataFrame(
        {
            "TCR-IR_VDJ_1_junction_aa": ["CASSA", "CASST"],
            "TCR-IR_VDJ_1_v_call": ["TRBV1", "TRBV2"],
        },
        index=["not-c1", "not-c2"],
    )
    adata = DummyAdata(obs_names=["c1", "c2"], uns={"TCR_IR": wells})

    with pytest.raises(ValueError, match="refusing ambiguous row alignment"):
        extract_wells_tcr_ir_features(adata)


def test_extract_wells_tcr_ir_rejects_length_mismatch():
    from scrfu.extract import extract_wells_tcr_ir_features

    wells = {
        "TCR-IR_VDJ_1_junction_aa": ["CASSA"],
        "TCR-IR_VDJ_1_v_call": ["TRBV1"],
    }
    adata = DummyAdata(obs_names=["c1", "c2"], uns={"TCR_IR": wells})

    with pytest.raises(ValueError, match="1 table rows versus 2 observations"):
        extract_wells_tcr_ir_features(adata)


def test_extract_wells_tcr_ir_rejects_duplicate_explicit_cells():
    from scrfu.extract import extract_wells_tcr_ir_features

    wells = pd.DataFrame(
        {
            "cell_id": ["c1", "c1"],
            "TCR-IR_VDJ_1_junction_aa": ["CASSA", "CASST"],
            "TCR-IR_VDJ_1_v_call": ["TRBV1", "TRBV2"],
        }
    )
    adata = DummyAdata(obs_names=["c1"], uns={"TCR_IR": wells})

    with pytest.raises(ValueError, match="duplicate cell identifiers"):
        extract_wells_tcr_ir_features(adata)


def test_extract_wells_tcr_ir_reports_missing_flattened_columns():
    from scrfu.extract import extract_wells_tcr_ir_features

    wells = pd.DataFrame({"TCR-IR_VDJ_1_junction_aa": ["CASSA"]}, index=["c1"])
    adata = DummyAdata(obs_names=["c1"], uns={"TCR_IR": wells})

    with pytest.raises(ValueError, match="VDJ_1 V call"):
        extract_wells_tcr_ir_features(adata)
