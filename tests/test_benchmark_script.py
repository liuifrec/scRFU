from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = "examples/benchmark_scirpy_dataset.py"


def _write_synthetic_h5ad(path: Path) -> None:
    obs_names = pd.Index(["c1", "c2", "c3"])
    obs = pd.DataFrame(
        {
            "sample": ["s1", "s1", "s2"],
            "cell_type": ["CD8_T", "CD4_T", "CD8_T"],
            "rfu_label": ["RFU1", "RFU2", "RFU1"],
            "rfu_score": [0.8, 0.7, 0.9],
        },
        index=obs_names,
    )
    adata = ad.AnnData(obs=obs)
    adata.obsm["airr"] = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3"],
            "chain": ["TRB", "TRB", "TRB"],
            "cdr3aa": ["CASSA", "CASST", "CASSG"],
            "v_call": ["TRBV1", "TRBV2", "TRBV1"],
            "productive": [True, True, True],
        },
        index=obs_names,
    )
    adata.write_h5ad(path)


def test_benchmark_help_works() -> None:
    result = subprocess.run(
        [sys.executable, SCRIPT, "--help"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert "usage:" in result.stdout


def test_benchmark_dry_run_missing_input_fails_clearly(tmp_path: Path) -> None:
    result = subprocess.run(
        [
            sys.executable,
            SCRIPT,
            "--input",
            str(tmp_path / "missing.h5ad"),
            "--rfu-dir",
            str(tmp_path / "RFU"),
            "--groupby",
            "sample",
            "--outdir",
            str(tmp_path / "out"),
            "--dry-run",
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "input file not found" in result.stderr


def test_benchmark_dry_run_with_synthetic_h5ad_succeeds(tmp_path: Path) -> None:
    input_path = tmp_path / "input.h5ad"
    outdir = tmp_path / "out"
    _write_synthetic_h5ad(input_path)

    result = subprocess.run(
        [
            sys.executable,
            SCRIPT,
            "--input",
            str(input_path),
            "--rfu-dir",
            str(tmp_path / "RFU"),
            "--groupby",
            "sample",
            "--outdir",
            str(outdir),
            "--dry-run",
            "--max-cells",
            "2",
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert (outdir / "validate_airr.tsv").exists()
    manifest = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["rfu_run"] is False
    assert manifest["n_cells"] == 2


def test_benchmark_skip_rfu_with_synthetic_h5ad_writes_outputs(tmp_path: Path) -> None:
    input_path = tmp_path / "input.h5ad"
    outdir = tmp_path / "out"
    _write_synthetic_h5ad(input_path)

    result = subprocess.run(
        [
            sys.executable,
            SCRIPT,
            "--input",
            str(input_path),
            "--rfu-dir",
            str(tmp_path / "RFU"),
            "--groupby",
            "sample",
            "--cell-type-key",
            "cell_type",
            "--outdir",
            str(outdir),
            "--skip-rfu",
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    for filename in [
        "validate_airr.tsv",
        "rfu_summary_global.tsv",
        "rfu_summary_by_group.tsv",
        "rfu_matrix_by_group.tsv",
        "rfu_bar_by_group.png",
        "rfu_heatmap_by_group.png",
        "rfu_score_hist.png",
        "rfu_summary_by_cell_type.tsv",
        "rfu_matrix_by_cell_type.tsv",
        "rfu_heatmap_by_cell_type.png",
        "run_manifest.json",
    ]:
        assert (outdir / filename).exists(), filename

    manifest = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["rfu_run"] is False
    assert manifest["groupby"] == "sample"
