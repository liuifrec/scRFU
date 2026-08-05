from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "examples" / "synthetic_scirpy_demo.py"


def _subprocess_output(result: subprocess.CompletedProcess[str]) -> str:
    return f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"


def test_synthetic_scirpy_demo_runs_from_unrelated_working_directory(tmp_path: Path) -> None:
    demo_outdir = tmp_path / "demo"
    env = os.environ.copy()
    env["SCRFU_DEMO_OUTDIR"] = str(demo_outdir)

    result = subprocess.run(
        [sys.executable, str(SCRIPT)],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, _subprocess_output(result)
    assert (demo_outdir / "rfu_matrix.tsv").exists()
    assert (demo_outdir / "rfu_bar.png").exists()
    assert (demo_outdir / "rfu_heatmap.png").exists()
    assert (demo_outdir / "rfu_score_hist.png").exists()
