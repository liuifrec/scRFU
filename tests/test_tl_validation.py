from __future__ import annotations

import pandas as pd
import pytest

from scrfu.tl import call_rfu
from scrfu.tl_rfu_repo import call_rfu_repo


class DummyAdata:
    def __init__(self) -> None:
        self.obsm = {
            "airr": pd.DataFrame(
                {
                    "cell_id": ["c1"],
                    "chain": ["TRB"],
                    "cdr3aa": ["CASSA"],
                    "v_call": ["TRBV1"],
                    "productive": [True],
                }
            )
        }
        self.obs_names = ["c1"]
        self.obs = pd.DataFrame(index=pd.Index(self.obs_names))
        self.uns = {}


def test_call_rfu_rejects_unknown_backend():
    with pytest.raises(ValueError, match="Unknown backend"):
        call_rfu(object(), backend="unknown")


def test_call_rfu_requires_rfu_dir_for_rfu_repo_backend():
    with pytest.raises(ValueError, match="requires rfu_dir"):
        call_rfu(object(), backend="rfu_repo")


def test_call_rfu_rejects_string_extra_r_args():
    with pytest.raises(TypeError, match="sequence of strings"):
        call_rfu(object(), backend="rfu_repo", rfu_dir="~/ext/RFU", extra_r_args="--vanilla")


def test_call_rfu_repo_rejects_string_extra_r_args():
    with pytest.raises(TypeError, match="sequence of strings"):
        call_rfu_repo(DummyAdata(), rfu_dir="~/ext/RFU", extra_r_args="--vanilla")


def test_call_rfu_repo_surfaces_missing_rfu_repo_files(tmp_path):
    with pytest.raises(FileNotFoundError, match="RFU_DIR is missing required files"):
        call_rfu_repo(DummyAdata(), rfu_dir=tmp_path)
