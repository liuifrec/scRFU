import pandas as pd
import pytest


class DummyAdata:
    def __init__(self, obs_names):
        self.obs_names = obs_names
        self.obs = pd.DataFrame(index=pd.Index(obs_names))
        self.uns = {}


def test_attach_rfu_results():
    from scrfu.attach import attach_rfu_results

    adata = DummyAdata(["c1", "c2", "c3"])

    features = pd.DataFrame(
        {"cell_id": ["c1", "c3"], "cdr3aa": ["CASSA", "CASST"], "trbv": ["TRBV1", "TRBV2"]}
    )
    rfu = pd.DataFrame(
        {"cell_id": ["c1", "c3"], "rfu_label": ["RFU3386", "RFU858"], "rfu_score": [0.1, 0.2]}
    )

    attach_rfu_results(adata, features, rfu)

    assert "rfu_label" in adata.obs.columns
    assert adata.obs.loc["c1", "rfu_label"] == "RFU3386"
    assert pd.isna(adata.obs.loc["c2", "rfu_label"])
    assert adata.uns["scrfu"]["package_version"]


def test_attach_rfu_results_accepts_alternate_columns_and_records_provenance():
    from scrfu.attach import attach_rfu_results

    adata = DummyAdata(["c1", "c2"])
    features = pd.DataFrame({"cell_id": ["c1"], "cdr3_aa": ["CASSA"], "v_gene": ["TRBV1"]})
    rfu = pd.DataFrame({"cell_id": ["c1"], "label": ["RFU7"], "score": ["0.25"]})

    attach_rfu_results(adata, features, rfu, provenance={"source": "unit-test"})

    assert adata.obs.loc["c1", "trb_cdr3aa"] == "CASSA"
    assert adata.obs.loc["c1", "trbv"] == "TRBV1"
    assert adata.obs.loc["c1", "rfu_label"] == "RFU7"
    assert adata.obs.loc["c1", "rfu_score"] == pytest.approx(0.25)
    assert adata.uns["scrfu"]["rfu"]["source"] == "unit-test"


def test_attach_rfu_results_handles_empty_features():
    from scrfu.attach import attach_rfu_results

    adata = DummyAdata(["c1", "c2"])
    features = pd.DataFrame(columns=["cell_id", "cdr3aa", "trbv"])

    attach_rfu_results(adata, features, pd.DataFrame(), provenance={"note": "empty"})

    assert list(adata.obs.columns) == ["trb_cdr3aa", "trbv", "rfu_label", "rfu_score"]
    assert adata.obs["rfu_score"].isna().all()
    assert adata.uns["scrfu"]["rfu"]["note"] == "empty"
