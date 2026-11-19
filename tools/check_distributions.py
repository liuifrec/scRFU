"""CI gate: separate clean wheel/sdist installs, checked outside the repository."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    for pattern in ["*.whl", "*.tar.gz"]:
        artifacts = list((root / "dist").glob(pattern))
        if len(artifacts) != 1:
            raise RuntimeError(f"Expected exactly one {pattern} artifact, found {artifacts}")
        artifact = artifacts[0]
        if artifact.suffix == ".whl":
            with zipfile.ZipFile(artifact) as archive:
                members = archive.namelist()
            assert "scrfu/resources/run_rfu_repo.R" in members
            assert "scrfu/data/tutorial_receptors.tsv" in members
        else:
            with tarfile.open(artifact) as archive:
                members = [name.partition("/")[2] for name in archive.getnames()]
            assert "src/scrfu/resources/run_rfu_repo.R" in members
            assert "tools/run_installed_tests.py" in members
        forbidden = ("manuscript/", "notebooks/", ".scrfu/", "results/")
        assert not any(
            name.startswith(forbidden) or name.endswith((".Rdata", ".h5ad", ".ipynb"))
            for name in members
        )
        with tempfile.TemporaryDirectory(prefix="scrfu-installed-") as temporary:
            stage = Path(temporary)
            subprocess.run([sys.executable, "-m", "venv", str(stage / "venv")], check=True)
            python = stage / "venv" / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
            subprocess.run(
                [str(python), "-m", "pip", "install", str(artifacts[0])], check=True, env=env
            )
            outside = stage / "outside"
            outside.mkdir()
            shutil.copy2(root / "tools/verify_installed.py", outside / "verify_installed.py")
            subprocess.run(
                [
                    str(python),
                    "verify_installed.py",
                    "--source-root",
                    str(root),
                    "--require-core-only",
                ],
                cwd=outside,
                env=env,
                check=True,
            )
            subprocess.run([str(python), "-m", "pip", "check"], check=True, env=env)
            # Upgrade only the optional plotting extra after proving core isolation.
            subprocess.run(
                [str(python), "-m", "pip", "install", f"scrfu[plotting] @ {artifact.as_uri()}"],
                check=True,
                env=env,
            )
            shutil.copy2(root / "examples/vdj_plotting.py", outside / "vdj_plotting.py")
            subprocess.run(
                [str(python), "vdj_plotting.py", "--outdir", "gallery"],
                cwd=outside,
                env=env,
                check=True,
            )
            subprocess.run(
                [
                    str(python),
                    "-c",
                    "import importlib.util; assert all(importlib.util.find_spec(x) is None for x in ['scirpy', 'scanpy', 'mudata', 'networkx', 'logomaker'])",
                ],
                cwd=outside,
                env=env,
                check=True,
            )
            subprocess.run([str(python), "-m", "pip", "check"], check=True, env=env)
            print(
                f"PASS: {artifact.name}: core-only then plotting-extra gallery outside checkout",
                flush=True,
            )


if __name__ == "__main__":
    main()
