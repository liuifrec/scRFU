from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from scrfu.completed_run import validate_completed_rfu_run


def _write_run(path: Path) -> Path:
    path.mkdir()
    receptors = pd.DataFrame(
        {
            "input_row_id": ["r0", "r1", "r2"],
            "cell_id": ["c0", "c1", "c2"],
            "chain": ["TRB", "TRB", "TRB"],
            "cdr3aa": ["CASSA", "CASSA", "BAD"],
            "v_call": ["TRBV1", "TRBV2", "TRBV3"],
        }
    )
    mapping = receptors[["input_row_id", "cell_id", "cdr3aa"]].copy()
    mapping["unique_sequence_id"] = ["sequence_00000000", "sequence_00000000", pd.NA]
    mapping["trbv"] = receptors["v_call"]
    mapping["eligibility_status"] = ["eligible", "eligible", "ineligible_cdr3_not_starting_c"]
    sequences = pd.DataFrame(
        {
            "unique_sequence_id": ["sequence_00000000"],
            "cdr3aa": ["CASSA"],
            "chain": ["TRB"],
            "v_call": [pd.NA],
            "rfu_id": [1],
            "rfu_label": ["RFU1"],
            "rfu_score": [0.8],
            "pass_thr": [True],
        }
    )
    rows = receptors.copy()
    rows["eligibility_status"] = mapping["eligibility_status"]
    rows["unique_sequence_id"] = mapping["unique_sequence_id"]
    rows["rfu_id"] = [1, 1, pd.NA]
    rows["rfu_label"] = ["RFU1", "RFU1", pd.NA]
    rows["rfu_score"] = [0.8, 0.8, pd.NA]
    rows["pass_thr"] = [True, True, pd.NA]
    rows["assignment_status"] = [
        "nearest_threshold_qualified",
        "nearest_threshold_qualified",
        "ineligible_sequence",
    ]
    receptors.to_csv(path / "receptors.tsv.gz", sep="\t", index=False)
    mapping.to_csv(path / "unique_sequence_map.tsv.gz", sep="\t", index=False)
    sequences.to_csv(path / "rfu_results_per_sequence.tsv.gz", sep="\t", index=False)
    rows.to_csv(path / "rfu_results_per_row.tsv.gz", sep="\t", index=False)
    (path / "run_manifest.json").write_text(
        json.dumps(
            {
                "original_row_count": 3,
                "rfu_threshold": 0.6,
                "chunk_count": 0,
                "run_manifest_path": None,
            }
        )
    )
    return path


def test_completed_run_validator_accepts_v_heterogeneity(tmp_path: Path) -> None:
    run = _write_run(tmp_path / "run")
    report = validate_completed_rfu_run(run)
    assert report["status"] == "valid"
    assert report["row_count"] == 3
    assert report["unique_sequence_count"] == 1
    assert report["threshold_pass_count"] == 2


def test_completed_run_validator_detects_reconstruction_order_change(tmp_path: Path) -> None:
    run = _write_run(tmp_path / "run")
    rows = pd.read_csv(run / "rfu_results_per_row.tsv.gz", sep="\t")
    rows.iloc[::-1].to_csv(run / "rfu_results_per_row.tsv.gz", sep="\t", index=False)
    report = validate_completed_rfu_run(run, strict=False)
    assert report["status"] == "invalid"
    assert any("preserve input order" in error for error in report["errors"])


def test_completed_run_validator_detects_bad_threshold_status(tmp_path: Path) -> None:
    run = _write_run(tmp_path / "run")
    rows = pd.read_csv(run / "rfu_results_per_row.tsv.gz", sep="\t")
    rows.loc[0, "assignment_status"] = "nearest_below_threshold"
    rows.to_csv(run / "rfu_results_per_row.tsv.gz", sep="\t", index=False)
    report = validate_completed_rfu_run(run, strict=False)
    assert any("assignment_status" in error for error in report["errors"])
