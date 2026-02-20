import pandas as pd


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