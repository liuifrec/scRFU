#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import resource
import time
from pathlib import Path

import pandas as pd


def _read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t" if ".tsv" in path.name else ",")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def compare(
    scrfu_path: Path,
    original_path: Path,
    *,
    id_col: str,
    score_tolerance: float,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    left, right = _read(scrfu_path), _read(original_path)
    required = [id_col, "rfu_id", "rfu_score", "pass_thr"]
    for label, frame in (("scRFU", left), ("original RFU", right)):
        missing = [column for column in required if column not in frame]
        if missing:
            raise ValueError(f"{label} table is missing columns: {missing}")
        if frame[id_col].isna().any() or frame[id_col].duplicated().any():
            raise ValueError(f"{label} {id_col!r} values must be unique and non-missing.")
    if set(left[id_col]) != set(right[id_col]):
        only_left = sorted(set(left[id_col]).difference(right[id_col]), key=str)
        only_right = sorted(set(right[id_col]).difference(left[id_col]), key=str)
        raise ValueError(
            f"Identifier sets differ; scRFU-only={only_left[:5]}, original-only={only_right[:5]}"
        )
    right = right.set_index(id_col).reindex(left[id_col]).reset_index()
    output = pd.DataFrame({id_col: left[id_col]})
    for column in ("rfu_id", "rfu_score", "pass_thr"):
        output[f"scrfu_{column}"] = left[column].to_numpy()
        output[f"original_{column}"] = right[column].to_numpy()
    output["rfu_id_match"] = (
        output["scrfu_rfu_id"]
        .astype("string")
        .fillna("<NA>")
        .eq(output["original_rfu_id"].astype("string").fillna("<NA>"))
    )
    left_score = pd.to_numeric(output["scrfu_rfu_score"], errors="coerce")
    right_score = pd.to_numeric(output["original_rfu_score"], errors="coerce")
    output["score_absolute_error"] = (left_score - right_score).abs()
    output["rfu_score_match"] = output["score_absolute_error"].le(score_tolerance) | (
        left_score.isna() & right_score.isna()
    )
    output["pass_thr_match"] = (
        output["scrfu_pass_thr"].astype("boolean").eq(output["original_pass_thr"].astype("boolean"))
    )
    output["row_match"] = output[["rfu_id_match", "rfu_score_match", "pass_thr_match"]].all(axis=1)
    mismatches = output.loc[~output["row_match"]].copy()
    summary: dict[str, object] = {
        "schema_version": "1",
        "row_count": len(output),
        "mismatch_count": len(mismatches),
        "rfu_id_mismatch_count": int((~output["rfu_id_match"]).sum()),
        "rfu_score_mismatch_count": int((~output["rfu_score_match"]).sum()),
        "threshold_mismatch_count": int((~output["pass_thr_match"]).sum()),
        "maximum_score_absolute_error": float(output["score_absolute_error"].max())
        if output["score_absolute_error"].notna().any()
        else None,
        "score_tolerance": score_tolerance,
        "passed": not len(mismatches),
        "scrfu_input_sha256": _sha256(scrfu_path),
        "original_input_sha256": _sha256(original_path),
    }
    return output, mismatches, summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare precomputed scRFU and frozen original-RFU row outputs."
    )
    parser.add_argument("--scrfu-results", type=Path, required=True)
    parser.add_argument("--original-results", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--id-col", default="input_row_id")
    parser.add_argument("--score-tolerance", type=float, default=1e-12)
    args = parser.parse_args()
    if args.score_tolerance < 0:
        parser.error("--score-tolerance must be non-negative")
    started = time.perf_counter()
    comparison, mismatches, summary = compare(
        args.scrfu_results,
        args.original_results,
        id_col=args.id_col,
        score_tolerance=args.score_tolerance,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(args.outdir / "row_comparison.tsv.gz", sep="\t", index=False)
    mismatches.to_csv(args.outdir / "mismatches.tsv.gz", sep="\t", index=False)
    summary["elapsed_seconds"] = time.perf_counter() - started
    summary["peak_rss_kb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
