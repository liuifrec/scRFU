#!/usr/bin/env python3
"""Combine panel-level source manifests without changing frozen evidence."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from _figure_common import MANIFEST_COLUMNS


def main(args: argparse.Namespace) -> None:
    frames = [pd.read_csv(path, sep="\t", dtype=str) for path in args.manifest]
    combined = pd.concat(frames, ignore_index=True)
    if list(combined.columns) != MANIFEST_COLUMNS:
        raise ValueError("A panel manifest does not match the frozen manifest schema.")
    if combined[["figure", "panel"]].duplicated().any():
        raise ValueError("Duplicate figure/panel provenance entries detected.")
    if combined[MANIFEST_COLUMNS].isna().any().any():
        raise ValueError("Every generated panel must have complete provenance.")
    combined.to_csv(args.output, sep="\t", index=False)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--manifest", type=Path, action="append", required=True)
    result.add_argument("--output", type=Path, required=True)
    return result


if __name__ == "__main__":
    main(parser().parse_args())
