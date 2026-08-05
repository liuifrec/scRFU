from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_DIR = REPO_ROOT / "examples"


def _subprocess_output(result: subprocess.CompletedProcess[str]) -> str:
    return f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"


@pytest.mark.parametrize(
    "script",
    [
        "real_scirpy_workflow.py",
        "radiation_pbmc_workflow.py",
    ],
)
def test_workflow_script_help_returns_zero(script: str, tmp_path: Path) -> None:
    result = subprocess.run(
        [sys.executable, str(EXAMPLES_DIR / script), "--help"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, _subprocess_output(result)
    assert "usage:" in result.stdout, _subprocess_output(result)


@pytest.mark.parametrize(
    "script",
    [
        "real_scirpy_workflow.py",
        "radiation_pbmc_workflow.py",
    ],
)
def test_workflow_scripts_fail_cleanly_when_input_missing(script: str, tmp_path: Path) -> None:
    missing_input = tmp_path / "missing.h5ad"
    result = subprocess.run(
        [
            sys.executable,
            str(EXAMPLES_DIR / script),
            "--input",
            str(missing_input),
            "--rfu-dir",
            str(tmp_path),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0, _subprocess_output(result)
    assert "input file not found" in result.stderr, _subprocess_output(result)


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
        "real_scirpy_workflow.py",
        "radiation_pbmc_workflow.py",
    ],
)
def test_workflow_examples_exist(path: str) -> None:
    assert (EXAMPLES_DIR / path).is_file()
