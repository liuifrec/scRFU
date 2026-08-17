from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

import scrfu

ROOT = Path(__file__).parents[1]


def _reference() -> scrfu.tl.FrozenRFUReference:
    return scrfu.tl.FrozenRFUReference.create(
        rfu_r_sha256="a" * 64,
        embedding_sha256="b" * 64,
        centroid_atlas_sha256="c" * 64,
        threshold=0.6,
        eligibility_rule="starts_with_C",
        assignment_mode="standard",
        receptor_chain="TRB",
        receptor_model="original_rfu_trb",
        reference_label="synthetic",
        creation_timestamp="2026-01-01T00:00:00+00:00",
    )


def _results() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "input_row_id": ["r1", "r2"],
            "cell_id": ["c1", "c2"],
            "sample": ["s1", "s2"],
            "cdr3aa": ["CASSF", "CASRF"],
            "rfu_id": [1, 2],
            "rfu_label": ["RFU1", "RFU2"],
            "rfu_score": [0.8, 0.7],
            "pass_thr": [True, True],
            "reference_coverage_status": ["covered", "covered"],
        }
    )


def test_original_rfu_parity_workflow(tmp_path: Path) -> None:
    left = tmp_path / "left.tsv"
    right = tmp_path / "right.tsv"
    _results().to_csv(left, sep="\t", index=False)
    changed = _results()
    changed.loc[1, "rfu_score"] += 0.1
    changed.to_csv(right, sep="\t", index=False)
    outdir = tmp_path / "parity"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "examples" / "original_rfu_parity.py"),
            "--scrfu-results",
            str(left),
            "--original-results",
            str(right),
            "--outdir",
            str(outdir),
        ],
        check=True,
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    summary = json.loads((outdir / "summary.json").read_text())
    assert summary["mismatch_count"] == 1
    assert not summary["passed"]
    assert (outdir / "mismatches.tsv.gz").is_file()


def test_cross_cohort_validation_workflow_and_help(tmp_path: Path) -> None:
    subprocess.run(
        [sys.executable, str(ROOT / "examples" / "cross_cohort_validation.py"), "--help"],
        check=True,
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    reference = _reference()
    cohorts = []
    for label, role in (("dev", "development"), ("hold", "held_out")):
        path = tmp_path / f"{label}.tsv"
        _results().to_csv(path, sep="\t", index=False)
        cohorts.append(
            {
                "label": label,
                "role": role,
                "results": str(path),
                "sample_key": "sample",
                "observed_reference_id": reference.immutable_reference_id,
            }
        )
    config = {
        "frozen_reference": reference.to_dict(),
        "frozen_parameters": {"threshold": 0.6},
        "data_hashes": {"dev": "d" * 64, "hold": "e" * 64},
        "evaluation_metrics": ["coverage"],
        "cohorts": cohorts,
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))
    outdir = tmp_path / "transfer"
    completed = subprocess.run(
        [
            sys.executable,
            str(ROOT / "examples" / "cross_cohort_validation.py"),
            "--config",
            str(config_path),
            "--outdir",
            str(outdir),
        ],
        check=True,
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert json.loads(completed.stdout)["held_out_cohort"] == "hold"
    assert (outdir / "heldout_manifest.json").is_file()
    summary = pd.read_csv(outdir / "cohort_summary.tsv", sep="\t")
    assert summary["cohort"].tolist() == ["dev", "hold"]
