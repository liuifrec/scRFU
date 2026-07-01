from __future__ import annotations

import pandas as pd
import pytest

from scrfu.io import export_rfu_matrix


class DummyAdata:
    def __init__(self, obs: pd.DataFrame) -> None:
        self.obs = obs


def test_export_rfu_matrix_tsv(tmp_path):
    adata = DummyAdata(
        pd.DataFrame(
            {
                "sample": ["s1", "s1", "s2"],
                "rfu_label": ["RFU1", "RFU2", "RFU1"],
            }
        )
    )

    path = tmp_path / "rfu.tsv"
    result = export_rfu_matrix(adata, groupby="sample", path=path)

    assert isinstance(result, pd.DataFrame)
    loaded = pd.read_csv(path, sep="\t", index_col=0)
    assert list(loaded.columns) == ["RFU1", "RFU2"]
    assert loaded.loc["s1"].sum() == pytest.approx(1.0)


def test_export_rfu_matrix_csv(tmp_path):
    adata = DummyAdata(
        pd.DataFrame(
            {
                "sample": ["s1", "s1", "s2"],
                "rfu_label": ["RFU1", "RFU2", "RFU1"],
            }
        )
    )

    path = tmp_path / "rfu.csv"
    result = export_rfu_matrix(adata, groupby="sample", path=path, sep=",")

    assert isinstance(result, pd.DataFrame)
    loaded = pd.read_csv(path, sep=",", index_col=0)
    assert list(loaded.columns) == ["RFU1", "RFU2"]


def test_export_rfu_matrix_creates_parent_dirs(tmp_path):
    adata = DummyAdata(pd.DataFrame({"sample": ["s1"], "rfu_label": ["RFU1"]}))

    path = tmp_path / "nested" / "rfu.tsv"
    export_rfu_matrix(adata, groupby="sample", path=path)

    assert path.exists()


def test_export_rfu_matrix_normalized_rows_sum_to_one(tmp_path):
    adata = DummyAdata(
        pd.DataFrame(
            {
                "sample": ["s1", "s1", "s2", "s2"],
                "rfu_label": ["RFU1", "RFU2", "RFU1", pd.NA],
            }
        )
    )

    result = export_rfu_matrix(adata, groupby="sample", path=tmp_path / "rfu.tsv", normalize=True)

    assert result.loc["s1"].sum() == pytest.approx(1.0)
    assert result.loc["s2"].sum() == pytest.approx(1.0)


def test_export_rfu_matrix_raw_counts(tmp_path):
    adata = DummyAdata(
        pd.DataFrame(
            {
                "sample": ["s1", "s1", "s1", "s2"],
                "rfu_label": ["RFU1", "RFU1", "RFU2", "RFU2"],
            }
        )
    )

    result = export_rfu_matrix(adata, groupby="sample", path=tmp_path / "rfu.tsv", normalize=False)

    assert result.loc["s1", "RFU1"] == 2
    assert result.loc["s1", "RFU2"] == 1
    assert result.loc["s2", "RFU2"] == 1


def test_export_rfu_matrix_missing_groupby_raises_value_error(tmp_path):
    adata = DummyAdata(pd.DataFrame({"rfu_label": ["RFU1"]}))

    with pytest.raises(ValueError, match="missing required columns"):
        export_rfu_matrix(adata, groupby="sample", path=tmp_path / "rfu.tsv")
