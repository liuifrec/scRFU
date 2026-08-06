from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import pytest

from examples.wells_atlas_workflow import extract_wells_features, main
from scrfu.rfu import RFURunResult

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "examples" / "wells_atlas_workflow.py"


def _synthetic_wells() -> ad.AnnData:
    obs_names = pd.Index(["c1", "c2", "c3"])
    adata = ad.AnnData(obs=pd.DataFrame(index=obs_names))
    adata.uns["TCR_IR"] = pd.DataFrame(
        {
            "TCR-IR_VDJ_1_locus": ["TRB", "TRB", "TRB"],
            "TCR-IR_VDJ_1_junction_aa": ["CASSA", "CASSA", "BAD"],
            "TCR-IR_VDJ_1_v_call": ["TRBV1", "TRBV2", "TRBV3"],
            "TCR-IR_VDJ_1_productive": [True, True, True],
            "TCR-IR_VDJ_2_locus": ["TRB", "TRA", "TRB"],
            "TCR-IR_VDJ_2_junction_aa": ["CASST", "CAVR", "nan"],
            "TCR-IR_VDJ_2_v_call": ["TRBV4", "TRAV1", "nan"],
            "TCR-IR_VDJ_2_productive": [True, True, False],
        },
        index=obs_names,
    )
    return adata


def _fake_rfu_dir(path: Path) -> Path:
    path.mkdir()
    (path / "RFU.R").write_text("AssignRFUs <- function(ff) {}\n", encoding="utf-8")
    (path / "trimerMDSfit_small.Rdata").write_bytes(b"trimer")
    (path / "km5000noMax.Rdata").write_bytes(b"centers")
    return path


def test_wells_workflow_help_runs_from_unrelated_directory(tmp_path: Path) -> None:
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--help"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    assert "--chunk-size" in result.stdout
    assert "--all-productive-chains" in result.stdout
    assert "--force-recompute" in result.stdout


def test_wells_all_productive_chain_extraction() -> None:
    features, qc = extract_wells_features(_synthetic_wells(), primary_chain=False)

    assert features[["cell_id", "cdr3aa", "chain_slot"]].to_dict("records") == [
        {"cell_id": "c1", "cdr3aa": "CASSA", "chain_slot": "VDJ_1"},
        {"cell_id": "c2", "cdr3aa": "CASSA", "chain_slot": "VDJ_1"},
        {"cell_id": "c3", "cdr3aa": "BAD", "chain_slot": "VDJ_1"},
        {"cell_id": "c1", "cdr3aa": "CASST", "chain_slot": "VDJ_2"},
    ]
    assert qc["extracted_trb_rows"] == 4
    assert qc["non_c_rows"] == 1


def test_wells_workflow_synthetic_smoke(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    input_path = tmp_path / "wells.h5ad"
    _synthetic_wells().write_h5ad(input_path)
    rfu_dir = _fake_rfu_dir(tmp_path / "RFU")
    outdir = tmp_path / "out"

    def fake_run(self, features: pd.DataFrame, **kwargs: object) -> RFURunResult:
        del self, kwargs
        rows = features.copy().reset_index(drop=True)
        rows["input_row_id"] = range(len(rows))
        rows["eligibility_status"] = [
            "eligible" if value.startswith("C") else "ineligible_cdr3_not_starting_c"
            for value in rows["cdr3aa"]
        ]
        rows["unique_sequence_id"] = pd.Series(
            ["sequence_00000000", "sequence_00000000", pd.NA], dtype="string"
        )
        rows["rfu_id"] = [10, 10, pd.NA]
        rows["rfu_label"] = ["RFU10", "RFU10", pd.NA]
        rows["rfu_score"] = [0.9, 0.9, pd.NA]
        rows["pass_thr"] = pd.Series([True, True, pd.NA], dtype="boolean")
        rows["rfu_status"] = [
            "assigned_threshold_pass",
            "assigned_threshold_pass",
            "ineligible_cdr3_not_starting_c",
        ]
        return RFURunResult(
            rows,
            "",
            "",
            0,
            metadata={
                "run_id": "synthetic-run",
                "chunking_enabled": True,
                "chunk_size": 2,
                "chunk_count": 1,
                "completed_chunk_count": 1,
                "reused_chunk_count": 0,
                "recomputed_chunk_count": 1,
                "invalidated_chunk_count": 0,
                "failed_chunk_count": 0,
                "total_elapsed_seconds": 0.01,
                "unique_query_count": 1,
                "original_row_count": 3,
                "reconstructed_output_row_count": 3,
                "run_manifest_path": None,
            },
        )

    monkeypatch.setattr("scrfu.backends.rfu_repo.RFURepoBackend.run", fake_run)

    returncode = main(
        [
            "--input",
            str(input_path),
            "--rfu-dir",
            str(rfu_dir),
            "--outdir",
            str(outdir),
            "--chunk-size",
            "2",
            "--primary-chain",
            "--write-annotated",
        ]
    )

    assert returncode == 0
    expected = {
        "extracted_trb.tsv.gz",
        "adapter_qc.json",
        "unique_sequence_map.tsv.gz",
        "rfu_results_per_sequence.tsv.gz",
        "rfu_results_per_cell.tsv.gz",
        "run_manifest.json",
        "annotated.h5ad",
    }
    assert expected.issubset(path.name for path in outdir.iterdir())
    manifest = json.loads((outdir / "run_manifest.json").read_text())
    assert manifest["backend_mode"] == "standard"
    assert manifest["run_id"] == "synthetic-run"
    per_cell = pd.read_csv(outdir / "rfu_results_per_cell.tsv.gz", sep="\t")
    assert per_cell["cell_id"].tolist() == ["c1", "c2", "c3"]
    sequence_map = pd.read_csv(outdir / "unique_sequence_map.tsv.gz", sep="\t")
    assert sequence_map["input_row_id"].tolist() == [0, 1, 2]
    assert sequence_map["cell_id"].tolist() == ["c1", "c2", "c3"]
    sequence_results = pd.read_csv(outdir / "rfu_results_per_sequence.tsv.gz", sep="\t")
    assert sequence_results["multiplicity"].tolist() == [2]
