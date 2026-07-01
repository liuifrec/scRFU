from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize(
    "script",
    [
        "examples/real_scirpy_workflow.py",
        "examples/radiation_pbmc_workflow.py",
    ],
)
def test_workflow_script_help_returns_zero(script: str) -> None:
    result = subprocess.run(
        [sys.executable, script, "--help"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert "usage:" in result.stdout


@pytest.mark.parametrize(
    "script",
    [
        "examples/real_scirpy_workflow.py",
        "examples/radiation_pbmc_workflow.py",
    ],
)
def test_workflow_scripts_fail_cleanly_when_input_missing(script: str, tmp_path: Path) -> None:
    missing_input = tmp_path / "missing.h5ad"
    result = subprocess.run(
        [
            sys.executable,
            script,
            "--input",
            str(missing_input),
            "--rfu-dir",
            str(tmp_path),
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "input file not found" in result.stderr


@pytest.mark.parametrize(
    "path",
    [
        "docs/manuscript_plan.md",
        "docs/roadmap.md",
        "docs/alfu_design_note.md",
        "docs/feature_table.md",
    ],
)
def test_manuscript_docs_exist(path: str) -> None:
    assert (REPO_ROOT / path).is_file()


@pytest.mark.parametrize(
    "path",
    [
        "examples/real_scirpy_workflow.py",
        "examples/radiation_pbmc_workflow.py",
    ],
)
def test_workflow_examples_exist(path: str) -> None:
    assert (REPO_ROOT / path).is_file()
