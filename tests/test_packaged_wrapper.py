"""Regression checks for installed and temporarily extracted wrapper lifetimes."""

from __future__ import annotations

import inspect
import pickle
import subprocess
import tempfile
from contextlib import contextmanager
from pathlib import Path

import pandas as pd
import pytest

from scrfu import _resources, tl
from scrfu.backends.rfu_repo import RFURepoBackend
from scrfu.cli import build_parser
from scrfu.doctor import doctor_report
from scrfu.tl_rfu_repo import call_rfu_repo, call_rfu_table


def _external(path):
    path.mkdir()
    for name, text in {
        "RFU.R": "AssignRFUs <- function(ff) {}",
        "trimerMDSfit_small.Rdata": "mock model",
        "km5000noMax.Rdata": "mock model",
    }.items():
        (path / name).write_text(text)
    return path


def _rows():
    return pd.DataFrame(
        {
            "input_row_id": ["r2", "r1", "r3"],
            "cell_id": ["b", "a", "c"],
            "chain": ["TRB"] * 3,
            "source_adapter": ["generic_airr_dataframe"] * 3,
            "cdr3aa": ["CASSLGQETQYF", "CASSIRSSYEQYF", "CASSLGQETQYF"],
            "v_call": ["TRBV1"] * 3,
        }
    )


@pytest.mark.parametrize("chunk_size", [None, 1])
def test_default_extraction_survives_subprocess_and_threads(tmp_path, monkeypatch, chunk_size):
    rfu_dir = _external(tmp_path / "external")
    monkeypatch.chdir(tmp_path)
    content = _resources._wrapper_resource().read_bytes()
    extracted = []
    invoked = []

    @contextmanager
    def temporary_resource(resource):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "wrapper.R"
            path.write_bytes(content)
            extracted.append(path)
            yield path

    def fake_run(cmd, **kwargs):
        # These are explicitly MOCK assignments. Existence is checked during use.
        wrapper = Path(cmd[1])
        assert wrapper.is_file() and wrapper.read_bytes() == content
        invoked.append(wrapper)
        queries = pd.read_csv(cmd[cmd.index("--in") + 1], sep="\t")
        pd.DataFrame(
            {
                "unique_sequence_id": queries.unique_sequence_id,
                "rfu_id": 1,
                "rfu_label": "RFU1",
                "rfu_score": 0.8,
                "pass_thr": True,
                "upstream_n_miss": 0,
            }
        ).to_csv(cmd[cmd.index("--out") + 1], sep="\t", index=False)
        return subprocess.CompletedProcess(cmd, 0, "mock", "")

    monkeypatch.setattr(_resources.resources, "as_file", temporary_resource)
    monkeypatch.setattr("scrfu.backends.rfu_repo.subprocess.run", fake_run)
    backend = pickle.loads(pickle.dumps(RFURepoBackend(rfu_dir=rfu_dir)))
    assert backend.wrapper_r_path is None
    result = call_rfu_table(
        _rows(),
        rfu_dir=rfu_dir,
        chunk_size=chunk_size,
        max_workers=2 if chunk_size else 1,
        executor="thread",
        workdir=tmp_path / "run",
    )
    assert result.per_row.input_row_id.tolist() == ["r2", "r1", "r3"]
    assert len(invoked) == (2 if chunk_size else 1)
    assert all(not path.exists() for path in extracted)
    assert result.provenance["wrapper_source"] == "package"


def test_explicit_override_and_doctor_agree(tmp_path):
    rfu_dir = _external(tmp_path / "external")
    explicit = tmp_path / "custom.R"
    explicit.write_text("# custom wrapper")
    backend = RFURepoBackend(rfu_dir=rfu_dir, wrapper_r_path=explicit)
    report = doctor_report(rfu_dir=rfu_dir, wrapper_r_path=explicit, verbose=True, environ={})
    assert backend.provenance_dict()["wrapper_r_path"] == report["wrapper"]["path"] == str(explicit)
    assert backend.provenance_dict()["wrapper_r_sha256"] == report["wrapper"]["sha256"]
    assert report["wrapper"]["source"] == "explicit"
    assert report["rfu_required_files_present"]
    assert not doctor_report(wrapper_r_path=tmp_path / "missing", environ={})["wrapper"][
        "available"
    ]


def test_every_public_default_uses_package_resolution():
    for function in [
        tl.call_rfu,
        tl.call_rfu_repo,
        tl.call_rfu_table,
        call_rfu_repo,
        call_rfu_table,
        tl.assign_rfu,
    ]:
        assert inspect.signature(function).parameters["wrapper_r_path"].default is None
    assert (
        build_parser().parse_args(["call-rfu", "in.h5ad", "-o", "out.h5ad"]).wrapper_r_path is None
    )
    assert doctor_report(environ={})["wrapper"]["source"] == "package"


def test_real_zip_resource_is_extracted_only_for_context(tmp_path, monkeypatch):
    import zipfile

    content = _resources._wrapper_resource().read_bytes()
    archive = tmp_path / "package.zip"
    with zipfile.ZipFile(archive, "w") as zipped:
        zipped.writestr("resources/run_rfu_repo.R", content)
    resource = zipfile.Path(archive, "resources/run_rfu_repo.R")
    monkeypatch.setattr(_resources, "_wrapper_resource", lambda: resource)
    with _resources.wrapper_path() as path:
        assert path.is_file() and path.read_bytes() == content
    assert not path.exists()
    assert _resources.wrapper_details()["path"] == "package:scrfu/resources/run_rfu_repo.R"
