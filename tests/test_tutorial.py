from __future__ import annotations

import json

import pandas as pd
import pytest

from examples.tutorial_end_to_end import run_tutorial


def test_synthetic_tutorial_end_to_end(tmp_path) -> None:
    manifest = run_tutorial(outdir=tmp_path)
    assert manifest["backend"] == "mock"
    assert all(value["synthetic"] for value in manifest["fixtures"].values())
    assert manifest["receptor_qc"]["status"] == "ok"
    assert not manifest["bcr_qc"]["tcr_rfu_assignment_permitted"]
    stored = json.loads((tmp_path / "run_manifest.json").read_text())
    assert stored["outputs"]["assigned_receptors.tsv"]["rows"] == 12
    reconstruction = pd.read_csv(tmp_path / "synthetic_vdjdb_row_summary.tsv", sep="\t")
    assert reconstruction["input_row_id"].tolist() == [f"r{index:02d}" for index in range(1, 13)]


def test_scverse_tutorial_when_optional_stack_is_installed(tmp_path) -> None:
    pytest.importorskip("scirpy")
    pytest.importorskip("mudata")
    from examples.tutorial_scverse_native import run

    output = run(tmp_path)
    assert output.is_file()


def test_external_backend_tutorial_preserves_sample_metadata(tmp_path, monkeypatch) -> None:
    from types import SimpleNamespace

    def mock_backend(receptors, **kwargs):
        assigned = receptors.copy()
        assigned["unique_sequence_id"] = assigned.cdr3aa
        assigned["rfu_id"] = 1
        assigned["rfu_label"] = "RFU1"
        assigned["rfu_score"] = 0.8
        assigned["pass_thr"] = True
        assigned["eligibility_status"] = "eligible"
        return SimpleNamespace(per_row=assigned, provenance={"test_backend": "mock"})

    monkeypatch.setattr("scrfu.tutorial.tl.call_rfu_table", mock_backend)
    manifest = run_tutorial(outdir=tmp_path, backend="rfu_repo", rfu_dir=tmp_path)
    assigned = pd.read_csv(tmp_path / "assigned_receptors.tsv", sep="\t")
    assert len(assigned) == 12 and assigned["sample"].nunique() == 6
    assert "sample_x" not in assigned and "sample_y" not in assigned
    assert manifest["outputs"]["rfu_pseudobulk_proportions.tsv"]["rows"] == 6
