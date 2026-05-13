from __future__ import annotations

from pathlib import Path

import pytest

from scrfu.backends.rfu_repo import RFURepoBackend, RFURepoPaths


def _write_required_rfu_files(rfu_dir: Path) -> None:
    for name in ("RFU.R", "trimerMDSfit_small.Rdata", "km5000noMax.Rdata"):
        (rfu_dir / name).write_text("placeholder", encoding="utf-8")


def test_rfu_repo_paths_from_dir_returns_expected_paths(tmp_path: Path):
    rfu_dir = tmp_path / "RFU"
    rfu_dir.mkdir()
    _write_required_rfu_files(rfu_dir)

    paths = RFURepoPaths.from_dir(rfu_dir)

    assert paths.rfu_dir == rfu_dir.resolve()
    assert paths.rfu_r == (rfu_dir / "RFU.R").resolve()
    assert paths.trimer_rdata == (rfu_dir / "trimerMDSfit_small.Rdata").resolve()
    assert paths.km5000_rdata == (rfu_dir / "km5000noMax.Rdata").resolve()


def test_rfu_repo_paths_from_dir_reports_missing_files(tmp_path: Path):
    rfu_dir = tmp_path / "RFU"
    rfu_dir.mkdir()
    (rfu_dir / "RFU.R").write_text("placeholder", encoding="utf-8")

    with pytest.raises(FileNotFoundError, match="trimerMDSfit_small.Rdata, km5000noMax.Rdata"):
        RFURepoPaths.from_dir(rfu_dir)


def test_backend_init_fails_when_rfu_dir_is_incomplete(tmp_path: Path):
    rfu_dir = tmp_path / "RFU"
    rfu_dir.mkdir()
    wrapper = tmp_path / "run_rfu_repo.R"
    wrapper.write_text("#!/usr/bin/env Rscript\n", encoding="utf-8")

    with pytest.raises(FileNotFoundError, match="RFU_DIR is missing required files"):
        RFURepoBackend(rfu_dir=rfu_dir, wrapper_r_path=wrapper)


def test_backend_init_fails_when_wrapper_script_is_missing(tmp_path: Path):
    rfu_dir = tmp_path / "RFU"
    rfu_dir.mkdir()
    _write_required_rfu_files(rfu_dir)

    with pytest.raises(FileNotFoundError, match="Wrapper R script not found"):
        RFURepoBackend(rfu_dir=rfu_dir, wrapper_r_path=tmp_path / "missing_wrapper.R")
