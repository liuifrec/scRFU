from __future__ import annotations

import os
import shutil
from pathlib import Path

import anndata as ad
import pandas as pd
import pytest

from scrfu.tl_rfu_repo import call_rfu_repo


def _integration_requirements() -> tuple[bool, str]:
    rscript = shutil.which("Rscript")
    if rscript is None:
        return False, "Rscript not available"

    rfu_dir_value = os.environ.get("RFU_DIR")
    if not rfu_dir_value:
        return False, "RFU_DIR is not set"

    rfu_dir = Path(rfu_dir_value).expanduser()
    required = ["RFU.R", "trimerMDSfit_small.Rdata", "km5000noMax.Rdata"]
    missing = [name for name in required if not (rfu_dir / name).exists()]
    if missing:
        return False, f"RFU_DIR is missing required files: {', '.join(missing)}"

    return True, ""


_READY, _REASON = _integration_requirements()


@pytest.mark.skipif(not _READY, reason=_REASON)
def test_call_rfu_repo_integration(tmp_path: Path):
    adata = ad.AnnData(obs=pd.DataFrame(index=pd.Index(["c1", "c2"])))
    adata.obsm["airr"] = pd.DataFrame(
        {
            "cell_id": ["c1", "c2"],
            "chain": ["TRB", "TRB"],
            "cdr3aa": ["CASSLGQETQYF", "CASSIRSSYEQYF"],
            "v_call": ["TRBV7-9", "TRBV6-5"],
            "productive": [True, True],
        }
    )

    result = call_rfu_repo(
        adata,
        rfu_dir=Path(os.environ["RFU_DIR"]).expanduser(),
        workdir=tmp_path / "work",
    )

    assert {"cell_id", "rfu_label", "rfu_score"}.issubset(result.columns)
    assert "trb_cdr3aa" in adata.obs.columns
    assert "trbv" in adata.obs.columns
    assert "rfu_label" in adata.obs.columns
    assert "rfu_score" in adata.obs.columns
    assert adata.uns["scrfu"]["package_version"]
    assert adata.uns["scrfu"]["rfu"]["backend"] == "rfu_repo"
