"""Seal existing runtime artifacts into a compact external evidence manifest.

This utility hashes files in place. It does not copy datasets or analysis
outputs into the source repository.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import scrfu


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _artifact(value: str) -> dict[str, object]:
    label, separator, raw_path = value.partition("=")
    if not separator or not label or not raw_path:
        raise ValueError("Artifacts must use LABEL=PATH syntax.")
    path = Path(raw_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Evidence artifact does not exist: {path}")
    return {
        "label": label,
        "runtime_only_path": str(path),
        "filename": path.name,
        "size_bytes": path.stat().st_size,
        "sha256": _sha256(path),
    }


def _git(command: list[str]) -> str | None:
    result = subprocess.run(["git", *command], capture_output=True, check=False, text=True)
    return result.stdout.strip() or None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-label", required=True)
    parser.add_argument("--analysis-type", required=True)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--input", action="append", default=[])
    parser.add_argument("--output", action="append", default=[])
    parser.add_argument("--source-manifest", action="append", default=[])
    parser.add_argument("--parameter", action="append", default=[])
    parser.add_argument("--random-seed", type=int)
    parser.add_argument("--runtime-seconds", type=float)
    parser.add_argument("--peak-rss-kb", type=int)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    parameters: dict[str, str] = {}
    for item in args.parameter:
        key, separator, value = item.partition("=")
        if not separator or not key:
            raise ValueError("Parameters must use KEY=VALUE syntax.")
        parameters[key] = value

    manifest = {
        "schema_version": 1,
        "dataset_label": args.dataset_label,
        "analysis_type": args.analysis_type,
        "public_source_accession": args.accession,
        "sealed_at": datetime.now(timezone.utc).isoformat(),
        "software": {
            "scrfu_version": scrfu.__version__,
            "git_commit": _git(["rev-parse", "HEAD"]),
            "git_dirty": bool(_git(["status", "--porcelain"])),
            "python_version": platform.python_version(),
            "platform": platform.platform(),
        },
        "parameters": parameters,
        "random_seed": args.random_seed,
        "runtime_seconds": args.runtime_seconds,
        "peak_rss_kb": args.peak_rss_kb,
        "inputs": [_artifact(value) for value in args.input],
        "source_manifests": [_artifact(value) for value in args.source_manifest],
        "outputs": [_artifact(value) for value in args.output],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"manifest": str(args.out), "sha256": _sha256(args.out)}, sort_keys=True))


if __name__ == "__main__":
    main()
