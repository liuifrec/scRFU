from __future__ import annotations

import os
import subprocess
import sys


def test_synthetic_scirpy_demo_runs_and_writes_outputs(tmp_path):
    demo_outdir = tmp_path / "demo"
    env = os.environ.copy()
    env["SCRFU_DEMO_OUTDIR"] = str(demo_outdir)

    result = subprocess.run(
        [sys.executable, "examples/synthetic_scirpy_demo.py"],
        cwd="/mnt/d/scRFU",
        env=env,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert (demo_outdir / "rfu_matrix.tsv").exists()
    assert (demo_outdir / "rfu_bar.png").exists()
    assert (demo_outdir / "rfu_heatmap.png").exists()
    assert (demo_outdir / "rfu_score_hist.png").exists()
