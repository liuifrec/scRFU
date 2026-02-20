import pandas as pd
import pytest


class DummyAdata:
    def __init__(self, airr_df: pd.DataFrame, obs_names):
        self.obsm = {"airr": airr_df}
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