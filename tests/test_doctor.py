from __future__ import annotations

from pathlib import Path

from scrfu.doctor import doctor_report


def _fake_rfu(path: Path) -> Path:
    path.mkdir()
    (path / "RFU.R").write_text("AssignRFUs <- function(ff) {}\n")
    (path / "trimerMDSfit_small.Rdata").write_bytes(b"trimer")
    (path / "km5000noMax.Rdata").write_bytes(b"centroids")
    return path


def test_doctor_unconfigured_is_read_only(tmp_path: Path) -> None:
    target = tmp_path / "not-created"
    report = doctor_report(output_dir=target, environ={})
    assert report["rfu_capability_mode"] == "unconfigured"
    assert report["rfu_artifact_hashes"] == {}
    assert report["output_writable"] is True
    assert not target.exists()


def test_doctor_valid_rfu_hashes_and_redacts(tmp_path: Path) -> None:
    rfu = _fake_rfu(tmp_path / "external-rfu")
    report = doctor_report(environ={"RFU_DIR": str(rfu)})
    assert report["rfu_required_files_present"] is True
    assert report["rfu_capability_mode"] == "standard"
    assert report["rfu_dir"] == "…/external-rfu"
    assert set(report["rfu_artifact_hashes"]) == {
        "RFU.R",
        "trimerMDSfit_small.Rdata",
        "km5000noMax.Rdata",
    }


def test_doctor_verbose_path_is_opt_in(tmp_path: Path) -> None:
    rfu = _fake_rfu(tmp_path / "external-rfu")
    report = doctor_report(verbose=True, environ={"RFU_DIR": str(rfu)})
    assert report["rfu_dir"] == str(rfu.resolve())
