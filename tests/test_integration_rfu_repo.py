from __future__ import annotations

import os
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from scrfu.backends.rfu_repo import RFUCapabilities, RFURepoBackend
from scrfu.tl_rfu_repo import call_rfu_repo

REPO_ROOT = Path(__file__).resolve().parents[1]
WRAPPER = REPO_ROOT / "r" / "run_rfu_repo.R"


def _checkout_requirements(environment_variable: str, capability: str) -> tuple[bool, str]:
    rscript = shutil.which("Rscript")
    if rscript is None:
        return False, "Rscript not available"

    rfu_dir_value = os.environ.get(environment_variable)
    if not rfu_dir_value:
        return False, f"{environment_variable} is not set"

    rfu_dir = Path(rfu_dir_value).expanduser()
    required = ["RFU.R", "trimerMDSfit_small.Rdata", "km5000noMax.Rdata"]
    missing = [name for name in required if not (rfu_dir / name).exists()]
    if missing:
        return False, f"{environment_variable} is missing required files: {', '.join(missing)}"

    capabilities = RFUCapabilities.from_rfu_r(rfu_dir / "RFU.R")
    if not getattr(capabilities, capability):
        function = {
            "assign_rfus": "AssignRFUs()",
            "assign_rfus_with_map": "AssignRFUs_with_map()",
        }[capability]
        return False, f"{environment_variable} does not provide {function}"

    return True, ""


_READY, _REASON = _checkout_requirements("RFU_DIR", "assign_rfus")
_MAP_READY, _MAP_REASON = _checkout_requirements("SCRFU_MAP_AWARE_RFU_DIR", "assign_rfus_with_map")


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
        },
        index=adata.obs_names,
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
    assert adata.uns["scrfu"]["rfu"]["backend_mode"] == "standard"
    assert adata.uns["scrfu"]["rfu"]["backend_resolution_source"] == "explicit_argument"


@pytest.mark.skipif(not _READY, reason=_REASON)
def test_unique_cdr3_scientific_equivalence_and_multiplicity_reconstruction(
    tmp_path: Path,
) -> None:
    features = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4", "c5", "c6", "c7"],
            "cdr3aa": [
                "CASSLGQETQYF",
                "CASSLGQETQYF",
                "CASSLGQETQYF",
                "CASSIRSSYEQYF",
                "CYYYYYYYYYYYF",
                "CSTQSTQSTQYF",
                "ASS_NOT_ELIGIBLE",
            ],
            "trbv": ["TRBV7-9", "TRBV1", "TRBV20-1", "TRBV6-5", "TRBV2", "TRBV3", "TRBV4"],
        }
    )
    backend = RFURepoBackend(
        rfu_dir=Path(os.environ["RFU_DIR"]),
        mode="standard",
        wrapper_r_path=WRAPPER,
    )

    repeated = backend.run(
        features,
        deduplicate=False,
        workdir=tmp_path / "repeated",
    )
    unique_features = (
        features[features["cdr3aa"].str.startswith("C")]
        .drop_duplicates("cdr3aa", keep="first")
        .reset_index(drop=True)
    )
    unique = backend.run(
        unique_features,
        deduplicate=False,
        workdir=tmp_path / "unique",
    )
    reconstructed = backend.run(
        features,
        deduplicate=True,
        workdir=tmp_path / "reconstructed",
    )

    eligible = repeated.df[repeated.df["eligibility_status"].eq("eligible")].reset_index(drop=True)
    rebuilt_eligible = reconstructed.df[
        reconstructed.df["eligibility_status"].eq("eligible")
    ].reset_index(drop=True)
    assert eligible["rfu_id"].tolist() == rebuilt_eligible["rfu_id"].tolist()
    assert eligible["rfu_label"].tolist() == rebuilt_eligible["rfu_label"].tolist()
    np.testing.assert_allclose(
        eligible["rfu_score"], rebuilt_eligible["rfu_score"], rtol=0.0, atol=1e-12
    )
    assert eligible["pass_thr"].tolist() == rebuilt_eligible["pass_thr"].tolist()

    unique_by_cdr3 = unique.df.set_index("cdr3aa")
    rebuilt_by_cdr3 = rebuilt_eligible.drop_duplicates("cdr3aa").set_index("cdr3aa")
    assert unique_by_cdr3["rfu_id"].to_dict() == rebuilt_by_cdr3["rfu_id"].to_dict()
    assert unique_by_cdr3["rfu_label"].to_dict() == rebuilt_by_cdr3["rfu_label"].to_dict()
    np.testing.assert_allclose(
        unique_by_cdr3.loc[rebuilt_by_cdr3.index, "rfu_score"],
        rebuilt_by_cdr3["rfu_score"],
        rtol=0.0,
        atol=1e-12,
    )
    assert unique_by_cdr3["pass_thr"].to_dict() == rebuilt_by_cdr3["pass_thr"].to_dict()

    assert (
        reconstructed.df["cdr3aa"].value_counts().to_dict()
        == features["cdr3aa"].value_counts().to_dict()
    )
    assert reconstructed.metadata["original_row_count"] == len(features)
    assert reconstructed.metadata["reconstructed_output_row_count"] == len(features)
    assert reconstructed.metadata["unique_query_count"] == len(unique_features)
    assert (
        repeated.metadata["upstream_threshold_miss_count"]
        == repeated.metadata["reconstructed_threshold_miss_count"]
    )
    assert (
        reconstructed.metadata["upstream_threshold_miss_count"]
        == unique.metadata["upstream_threshold_miss_count"]
    )
    assert (
        reconstructed.metadata["reconstructed_threshold_miss_count"]
        == repeated.metadata["upstream_threshold_miss_count"]
    )
    assert int((~eligible["pass_thr"].fillna(False)).sum()) >= 2


@pytest.mark.skipif(not _READY, reason=_REASON)
def test_official_standard_chunked_matches_unchunked(tmp_path: Path) -> None:
    features = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4", "c5"],
            "cdr3aa": [
                "CASSLGQETQYF",
                "CASSLGQETQYF",
                "CASSIRSSYEQYF",
                "CYYYYYYYYYYYF",
                "ASS_NOT_ELIGIBLE",
            ],
            "trbv": ["TRBV7-9", "TRBV1", "TRBV6-5", "TRBV2", "TRBV4"],
        }
    )
    backend = RFURepoBackend(
        rfu_dir=Path(os.environ["RFU_DIR"]),
        mode="standard",
        wrapper_r_path=WRAPPER,
    )

    unchunked = backend.run(features, workdir=tmp_path / "unchunked")
    chunked = backend.run(features, chunk_size=2, workdir=tmp_path / "chunked")

    columns = [
        "input_row_id",
        "unique_sequence_id",
        "cell_id",
        "rfu_id",
        "rfu_label",
        "rfu_score",
        "pass_thr",
        "eligibility_status",
        "rfu_status",
    ]
    pd.testing.assert_frame_equal(
        unchunked.df[columns],
        chunked.df[columns],
        check_exact=False,
        atol=1e-12,
        rtol=0.0,
    )


@pytest.mark.skipif(not _MAP_READY, reason=_MAP_REASON)
def test_map_aware_mode_does_not_trust_shifted_map_identifiers(tmp_path: Path) -> None:
    features = pd.DataFrame(
        {
            "cell_id": ["rejected", "later-1", "later-2"],
            "cdr3aa": ["ASS_REJECTED", "CASSLGQETQYF", "CASSIRSSYEQYF"],
            "trbv": ["TRBV1", "TRBV7-9", "TRBV6-5"],
        }
    )
    standard = RFURepoBackend(
        rfu_dir=Path(os.environ["SCRFU_MAP_AWARE_RFU_DIR"]),
        mode="standard",
        wrapper_r_path=WRAPPER,
    ).run(features, deduplicate=False, workdir=tmp_path / "standard")
    map_aware = RFURepoBackend(
        rfu_dir=Path(os.environ["SCRFU_MAP_AWARE_RFU_DIR"]),
        mode="map_aware",
        wrapper_r_path=WRAPPER,
    ).run(features, deduplicate=False, workdir=tmp_path / "map-aware")

    assert map_aware.df["input_row_id"].tolist() == [0, 1, 2]
    assert map_aware.df.loc[0, "rfu_status"] == "ineligible_cdr3_not_starting_c"
    assert pd.isna(map_aware.df.loc[0, "unique_sequence_id"])
    for row_id in (1, 2):
        assert map_aware.df.loc[row_id, "cell_id"] == features.loc[row_id, "cell_id"]
        assert map_aware.df.loc[row_id, "rfu_id"] == standard.df.loc[row_id, "rfu_id"]
        assert map_aware.df.loc[row_id, "rfu_label"] == standard.df.loc[row_id, "rfu_label"]
        assert map_aware.df.loc[row_id, "rfu_score"] == pytest.approx(
            standard.df.loc[row_id, "rfu_score"], abs=1e-12
        )
