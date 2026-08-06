from __future__ import annotations

import subprocess
from dataclasses import FrozenInstanceError
from pathlib import Path

import pandas as pd
import pytest

from scrfu._version import __version__
from scrfu.backends.rfu_repo import (
    RFUCapabilities,
    RFUCapabilityError,
    RFUConfigurationError,
    RFURepoBackend,
    RFURepoPaths,
)
from scrfu.rfu import file_sha256


def _write_fake_rfu_dir(
    path: Path,
    *,
    standard: bool = True,
    map_aware: bool = False,
    batch_maps: bool = False,
) -> Path:
    path.mkdir()
    definitions = []
    if standard:
        definitions.append("AssignRFUs <- function(ff) {}")
    if map_aware:
        definitions.append("AssignRFUs_with_map <- function(ff) {}")
    if batch_maps:
        definitions.append("RFUbatch_with_maps = function(DIR) {}")
    (path / "RFU.R").write_text("\n".join(definitions) + "\n", encoding="utf-8")
    (path / "trimerMDSfit_small.Rdata").write_bytes(b"synthetic trimer atlas")
    (path / "km5000noMax.Rdata").write_bytes(b"synthetic centroid atlas")
    return path


def _write_wrapper(path: Path) -> Path:
    path.write_text("#!/usr/bin/env Rscript\n", encoding="utf-8")
    return path


def _fake_successful_subprocess(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_run(cmd: list[str], **_: object) -> subprocess.CompletedProcess[str]:
        input_path = Path(cmd[cmd.index("--in") + 1])
        output_path = Path(cmd[cmd.index("--out") + 1])
        threshold = float(cmd[cmd.index("--threshold") + 1])
        queries = pd.read_csv(input_path, sep="\t")
        output = pd.DataFrame(
            {
                "unique_sequence_id": queries["unique_sequence_id"],
                "rfu_id": range(1, len(queries) + 1),
                "rfu_label": [f"RFU{i}" for i in range(1, len(queries) + 1)],
                "rfu_score": [0.8] * len(queries),
                "pass_thr": [threshold <= 0.8] * len(queries),
                "upstream_n_miss": [0] * len(queries),
            }
        )
        output.to_csv(output_path, sep="\t", index=False)
        return subprocess.CompletedProcess(cmd, 0, stdout="synthetic stdout", stderr="")

    monkeypatch.setattr("scrfu.backends.rfu_repo.subprocess.run", fake_run)


def test_explicit_rfu_dir_takes_precedence_over_environment(tmp_path: Path) -> None:
    explicit = _write_fake_rfu_dir(tmp_path / "explicit")
    environment = _write_fake_rfu_dir(tmp_path / "environment")

    paths = RFURepoPaths.resolve(explicit, environ={"RFU_DIR": str(environment)})

    assert paths.rfu_dir == explicit.resolve()
    assert paths.resolution_source == "explicit_argument"


def test_rfu_dir_falls_back_to_environment(tmp_path: Path) -> None:
    environment = _write_fake_rfu_dir(tmp_path / "environment")

    paths = RFURepoPaths.resolve(environ={"RFU_DIR": str(environment)})

    assert paths.rfu_dir == environment.resolve()
    assert paths.resolution_source == "environment_variable"


def test_missing_backend_configuration_has_actionable_error() -> None:
    with pytest.raises(RFUConfigurationError, match="rfu_dir='/path/to/RFU'.*RFU_DIR"):
        RFURepoPaths.resolve(environ={})


def test_rfu_repo_paths_validate_all_required_files(tmp_path: Path) -> None:
    rfu_dir = tmp_path / "RFU"
    rfu_dir.mkdir()
    (rfu_dir / "RFU.R").write_text("AssignRFUs <- function(ff) {}\n", encoding="utf-8")

    with pytest.raises(FileNotFoundError, match="trimerMDSfit_small.Rdata, km5000noMax.Rdata"):
        RFURepoPaths.from_dir(rfu_dir)


def test_capability_detection_is_static_exact_and_immutable(tmp_path: Path) -> None:
    rfu_r = tmp_path / "RFU.R"
    rfu_r.write_text(
        "# AssignRFUs_with_map <- function(ff) {}\n"
        "AssignRFUs <- function(ff) {}\n"
        "RFUbatch_with_maps = function(DIR) {}\n",
        encoding="utf-8",
    )

    capabilities = RFUCapabilities.from_rfu_r(rfu_r)

    assert capabilities.as_dict() == {
        "AssignRFUs": True,
        "AssignRFUs_with_map": False,
        "RFUbatch_with_maps": True,
    }
    with pytest.raises(FrozenInstanceError):
        capabilities.assign_rfus = False  # type: ignore[misc]

    without_batch = RFUCapabilities(True, False, False)
    with pytest.raises(RFUCapabilityError, match="RFUbatch_with_maps"):
        without_batch.require_batch_with_maps()


def test_backend_defaults_to_explicit_standard_mode(tmp_path: Path) -> None:
    rfu_dir = _write_fake_rfu_dir(tmp_path / "RFU", standard=True, map_aware=True)
    wrapper = _write_wrapper(tmp_path / "wrapper.R")

    backend = RFURepoBackend(rfu_dir=rfu_dir, wrapper_r_path=wrapper)

    assert backend.mode == "standard"
    assert backend.capabilities.assign_rfus
    assert backend.capabilities.assign_rfus_with_map


def test_backend_accepts_explicit_map_aware_mode(tmp_path: Path) -> None:
    rfu_dir = _write_fake_rfu_dir(tmp_path / "RFU", standard=True, map_aware=True)
    wrapper = _write_wrapper(tmp_path / "wrapper.R")

    backend = RFURepoBackend(rfu_dir=rfu_dir, mode="map_aware", wrapper_r_path=wrapper)

    assert backend.mode == "map_aware"


def test_backend_provenance_records_environment_resolution(tmp_path: Path) -> None:
    rfu_dir = _write_fake_rfu_dir(tmp_path / "RFU")
    wrapper = _write_wrapper(tmp_path / "wrapper.R")

    backend = RFURepoBackend(
        wrapper_r_path=wrapper,
        environ={"RFU_DIR": str(rfu_dir)},
    )

    assert backend.paths.resolution_source == "environment_variable"
    assert backend.provenance_dict()["backend_resolution_source"] == "environment_variable"


def test_map_aware_mode_requires_optional_function(tmp_path: Path) -> None:
    rfu_dir = _write_fake_rfu_dir(tmp_path / "RFU", standard=True)
    wrapper = _write_wrapper(tmp_path / "wrapper.R")

    with pytest.raises(RFUCapabilityError, match='Use mode="standard"'):
        RFURepoBackend(rfu_dir=rfu_dir, mode="map_aware", wrapper_r_path=wrapper)


def test_standard_mode_requires_public_assign_function(tmp_path: Path) -> None:
    rfu_dir = _write_fake_rfu_dir(tmp_path / "RFU", standard=False, map_aware=True)
    wrapper = _write_wrapper(tmp_path / "wrapper.R")

    with pytest.raises(RFUCapabilityError, match=r"AssignRFUs\(\)"):
        RFURepoBackend(rfu_dir=rfu_dir, wrapper_r_path=wrapper)


def test_backend_rejects_unknown_mode(tmp_path: Path) -> None:
    rfu_dir = _write_fake_rfu_dir(tmp_path / "RFU")
    wrapper = _write_wrapper(tmp_path / "wrapper.R")

    with pytest.raises(ValueError, match="Unknown RFU backend mode"):
        RFURepoBackend(rfu_dir=rfu_dir, mode="auto", wrapper_r_path=wrapper)


def test_backend_init_fails_when_wrapper_script_is_missing(tmp_path: Path) -> None:
    rfu_dir = _write_fake_rfu_dir(tmp_path / "RFU")

    with pytest.raises(FileNotFoundError, match="Wrapper R script not found"):
        RFURepoBackend(rfu_dir=rfu_dir, wrapper_r_path=tmp_path / "missing_wrapper.R")


def test_backend_provenance_records_resolution_capabilities_and_hashes(tmp_path: Path) -> None:
    rfu_dir = _write_fake_rfu_dir(tmp_path / "RFU", map_aware=True, batch_maps=True)
    wrapper = _write_wrapper(tmp_path / "wrapper.R")
    backend = RFURepoBackend(rfu_dir=rfu_dir, wrapper_r_path=wrapper)

    provenance = backend.provenance_dict()

    assert provenance["scrfu_version"] == __version__
    assert provenance["backend"] == "rfu_repo"
    assert provenance["backend_mode"] == "standard"
    assert provenance["backend_resolution_source"] == "explicit_argument"
    assert provenance["resolved_backend_path"] == str(rfu_dir.resolve())
    assert provenance["detected_capabilities"] == {
        "AssignRFUs": True,
        "AssignRFUs_with_map": True,
        "RFUbatch_with_maps": True,
    }
    for key in (
        "rfu_r_sha256",
        "trimer_rdata_sha256",
        "km5000_rdata_sha256",
        "wrapper_r_sha256",
    ):
        assert len(provenance[key]) == 64
    assert provenance["rfu_r_sha256"] == file_sha256(rfu_dir / "RFU.R")
    assert provenance["trimer_rdata_sha256"] == file_sha256(rfu_dir / "trimerMDSfit_small.Rdata")
    assert provenance["km5000_rdata_sha256"] == file_sha256(rfu_dir / "km5000noMax.Rdata")
    assert provenance["timestamp"].endswith("+00:00")


def test_standard_mapping_preserves_rows_order_and_duplicate_multiplicity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rfu_dir = _write_fake_rfu_dir(tmp_path / "RFU")
    wrapper = _write_wrapper(tmp_path / "wrapper.R")
    _fake_successful_subprocess(monkeypatch)
    backend = RFURepoBackend(rfu_dir=rfu_dir, wrapper_r_path=wrapper)
    features = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4", "c5"],
            "cdr3aa": ["CASSA", "NOT_C", "CASSA", "CASST", "CASSA"],
            "trbv": ["TRBV1", "TRBV2", "TRBV3", "TRBV4", "TRBV5"],
        }
    )

    run = backend.run(features, workdir=tmp_path / "work")

    assert run.df["input_row_id"].tolist() == [0, 1, 2, 3, 4]
    assert run.df["cell_id"].tolist() == features["cell_id"].tolist()
    assert len(run.df) == len(features)
    assert run.df.loc[[0, 2, 4], "unique_sequence_id"].nunique() == 1
    assert run.df.loc[[0, 2, 4], "rfu_id"].nunique() == 1
    assert run.df.loc[[0, 2, 4], "trbv"].tolist() == ["TRBV1", "TRBV3", "TRBV5"]
    assert run.df.loc[1, "eligibility_status"] == "ineligible_cdr3_not_starting_c"
    assert pd.isna(run.df.loc[1, "rfu_id"])
    expected_metadata = {
        "eligibility_rule": "^C",
        "deduplication_key": "cdr3aa",
        "original_row_count": 5,
        "eligible_row_count": 4,
        "unique_query_count": 2,
        "rfu_threshold": 0.6,
        "reconstructed_output_row_count": 5,
        "upstream_threshold_miss_count": 0,
        "reconstructed_threshold_miss_count": 0,
    }
    assert {key: run.metadata[key] for key in expected_metadata} == expected_metadata


def test_non_deduplicated_queries_preserve_one_query_per_eligible_row(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rfu_dir = _write_fake_rfu_dir(tmp_path / "RFU")
    wrapper = _write_wrapper(tmp_path / "wrapper.R")
    _fake_successful_subprocess(monkeypatch)
    backend = RFURepoBackend(rfu_dir=rfu_dir, wrapper_r_path=wrapper)
    features = pd.DataFrame(
        {
            "cell_id": ["c1", "c2"],
            "cdr3aa": ["CASSA", "CASSA"],
            "trbv": ["TRBV1", "TRBV2"],
        }
    )

    run = backend.run(features, deduplicate=False, workdir=tmp_path / "work")

    assert run.df["unique_sequence_id"].nunique() == 2
    assert run.metadata["deduplication_key"] is None
    assert run.metadata["unique_query_count"] == 2
