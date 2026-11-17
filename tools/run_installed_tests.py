"""Run software tests against an installed distribution, outside its checkout.

Copy test/example support files, never the Python implementation. No PYTHONPATH
or src-path injection is used. External R tests remain opt-in in standard CI.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import scrfu


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python-only", action="store_true")
    parser.add_argument("--junitxml", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    installed = Path(scrfu.__file__).resolve()
    assert installed.is_relative_to(Path(sys.prefix).resolve()), installed
    assert not installed.is_relative_to(root), installed
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    with tempfile.TemporaryDirectory(prefix="scrfu-test-support-") as temp:
        outside = Path(temp)
        for name in ["tests", "examples", "r", "docs", "benchmarks"]:
            shutil.copytree(
                root / name,
                outside / name,
                ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "_build"),
            )
        for name in ["pyproject.toml", "README.md", "CITATION.cff", "CHANGELOG.md", "LICENSE"]:
            shutil.copy2(root / name, outside / name)
        # Metadata integrity and the explicit source-checkout R launcher need
        # these two files, not an importable copy of the Python package.
        for name in ["src/scrfu/_version.py", "src/scrfu/resources/run_rfu_repo.R"]:
            target = outside / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(root / name, target)
        assert not (outside / "src/scrfu/__init__.py").exists()
        print(
            json.dumps(
                {
                    "package_file": str(installed),
                    "version": scrfu.__version__,
                    "python": sys.version,
                    "test_working_directory": str(outside),
                    "pythonpath": None,
                    "external_r_tests": "NOT RUN"
                    if args.python_only
                    else "configured by environment",
                },
                indent=2,
            ),
            flush=True,
        )
        command = [sys.executable, "-m", "pytest", "-ra"]
        if args.python_only:
            command.append("--ignore=tests/test_integration_rfu_repo.py")
        if args.junitxml:
            command.append(f"--junitxml={args.junitxml.resolve()}")
        subprocess.run(command, cwd=outside, env=env, check=True)


if __name__ == "__main__":
    main()
