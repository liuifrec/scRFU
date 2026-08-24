#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import resource
import subprocess
import sys
import time
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd


def _add_repo_src_to_path() -> None:
    src = Path(__file__).resolve().parents[1] / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))


def _read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t" if ".tsv" in path.name else ",")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _artifact_record(path: Path) -> dict[str, Any]:
    return {"sha256": _sha256(path), "size_bytes": path.stat().st_size}


def compare(
    scrfu_path: Path,
    original_path: Path,
    *,
    id_col: str,
    score_tolerance: float,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    left, right = _read(scrfu_path), _read(original_path)
    required = [id_col, "rfu_id", "rfu_label", "rfu_score", "pass_thr"]
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
    for column in ("rfu_id", "rfu_label", "rfu_score", "pass_thr"):
        output[f"scrfu_{column}"] = left[column].to_numpy()
        output[f"original_{column}"] = right[column].to_numpy()
    left_id = pd.to_numeric(output["scrfu_rfu_id"], errors="coerce")
    right_id = pd.to_numeric(output["original_rfu_id"], errors="coerce")
    output["rfu_id_match"] = left_id.eq(right_id) | (left_id.isna() & right_id.isna())
    output["rfu_label_match"] = (
        output["scrfu_rfu_label"]
        .astype("string")
        .fillna("<NA>")
        .eq(output["original_rfu_label"].astype("string").fillna("<NA>"))
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
    match_columns = ["rfu_id_match", "rfu_label_match", "rfu_score_match", "pass_thr_match"]
    output["row_match"] = output[match_columns].all(axis=1)
    mismatches = output.loc[~output["row_match"]].copy()
    summary: dict[str, object] = {
        "schema_version": "2",
        "row_count": len(output),
        "mismatch_count": len(mismatches),
        "exact_rfu_id_matches": int(output["rfu_id_match"].sum()),
        "rfu_id_mismatch_count": int((~output["rfu_id_match"]).sum()),
        "rfu_label_mismatch_count": int((~output["rfu_label_match"]).sum()),
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


def _fixture(seed: int) -> pd.DataFrame:
    rows = pd.DataFrame(
        {
            "input_row_id": [f"row_{index:02d}" for index in range(12)],
            "cell_id": [
                "unique",
                "duplicate_a",
                "duplicate_b",
                "different_v",
                "below_a",
                "below_b",
                "repeat_1",
                "repeat_2",
                "repeat_3",
                "second_unique",
                "ineligible_a",
                "ineligible_b",
            ],
            "cdr3aa": [
                "CASSLGQETQYF",
                "CASSIRSSYEQYF",
                "CASSIRSSYEQYF",
                "CASSIRSSYEQYF",
                "CYYYYYYYYYYYF",
                "CSTQSTQSTQYF",
                "CASSLGQETQYF",
                "CASSLGQETQYF",
                "CASSLGQETQYF",
                "CASSPGQGYEQYF",
                "ASS_NOT_ELIGIBLE",
                "GSQG_NOT_ELIGIBLE",
            ],
            "trbv": [
                "TRBV7-9",
                "TRBV6-5",
                "TRBV6-5",
                "TRBV20-1",
                "TRBV2",
                "TRBV3",
                "TRBV7-9",
                "TRBV7-9",
                "TRBV7-9",
                "TRBV5-1",
                "TRBV1",
                "TRBV4",
            ],
        }
    )
    return rows.sample(frac=1, random_state=seed).reset_index(drop=True)


def _git_value(checkout: Path, *args: str) -> str | None:
    result = subprocess.run(
        ["git", "-C", str(checkout), *args], capture_output=True, text=True, check=False
    )
    return result.stdout.strip() if result.returncode == 0 else None


def _run_official(
    eligible: pd.DataFrame,
    *,
    rfu_dir: Path,
    outdir: Path,
    threshold: float,
    rscript_bin: str,
) -> tuple[Path, float, int]:
    official_input = outdir / "official_input.tsv"
    official_output = outdir / "official_results.tsv"
    script_path = outdir / "run_official_assign_rfus.R"
    eligible.loc[:, ["input_row_id", "cdr3aa", "trbv"]].to_csv(
        official_input, sep="\t", index=False
    )
    script_path.write_text(
        """args <- commandArgs(trailingOnly=TRUE)
rfu_dir <- normalizePath(args[[1]])
input_path <- normalizePath(args[[2]])
output_path <- normalizePath(args[[3]], mustWork=FALSE)
threshold <- as.numeric(args[[4]])
source(file.path(rfu_dir, 'RFU.R'))
load(file.path(rfu_dir, 'trimerMDSfit_small.Rdata'))
load(file.path(rfu_dir, 'km5000noMax.Rdata'))
inp <- read.table(input_path, header=TRUE, sep='\t', quote='', stringsAsFactors=FALSE,
                  check.names=FALSE)
giana_path <- tempfile(fileext='.tsv')
giana <- data.frame(inp$cdr3aa, inp$trbv, 1L, 1L, inp$input_row_id)
write.table(giana, giana_path, sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
res <- AssignRFUs(giana_path, CL=km5000noMax, THR=threshold)
if (length(res$TCR) != nrow(inp) || length(res$COR) != nrow(inp)) {
  stop('official output length mismatch')
}
out <- data.frame(input_row_id=inp$input_row_id, rfu_id=as.integer(res$TCR),
                  rfu_label=paste0('RFU', as.integer(res$TCR)),
                  rfu_score=as.numeric(res$COR), pass_thr=as.numeric(res$COR) >= threshold)
write.table(out, output_path, sep='\t', quote=FALSE, row.names=FALSE)
""",
        encoding="utf-8",
    )
    started = time.perf_counter()
    completed = subprocess.run(
        [
            rscript_bin,
            str(script_path),
            str(rfu_dir),
            str(official_input),
            str(official_output),
            str(threshold),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    elapsed = time.perf_counter() - started
    official_peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if completed.returncode:
        raise RuntimeError(
            "Direct official RFU execution failed.\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    return official_output, elapsed, official_peak


def run_real_parity(
    *,
    rfu_dir: Path,
    outdir: Path,
    threshold: float,
    score_tolerance: float,
    seed: int,
    rscript_bin: str,
) -> dict[str, Any]:
    _add_repo_src_to_path()
    from scrfu.backends.rfu_repo import RFURepoBackend

    rfu_dir = rfu_dir.expanduser().resolve()
    outdir = outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    artifacts = {
        name: _artifact_record(rfu_dir / name)
        for name in ("RFU.R", "trimerMDSfit_small.Rdata", "km5000noMax.Rdata")
    }
    features = _fixture(seed)
    features.to_csv(outdir / "parity_fixture.tsv", sep="\t", index=False)
    eligible_mask = features["cdr3aa"].str.startswith("C", na=False)
    eligible = features.loc[eligible_mask].copy()

    backend = RFURepoBackend(
        rfu_dir=rfu_dir,
        mode="standard",
        wrapper_r_path=Path(__file__).resolve().parents[1] / "r" / "run_rfu_repo.R",
        rscript_bin=rscript_bin,
    )
    started = time.perf_counter()
    scrfu_run = backend.run(
        features,
        threshold=threshold,
        deduplicate=True,
        workdir=outdir / "scrfu_backend",
    )
    scrfu_elapsed = time.perf_counter() - started
    scrfu_results = outdir / "scrfu_results.tsv"
    scrfu_run.df.to_csv(scrfu_results, sep="\t", index=False)
    scrfu_eligible_results = outdir / "scrfu_eligible_results.tsv"
    scrfu_run.df.loc[scrfu_run.df["eligibility_status"].eq("eligible")].to_csv(
        scrfu_eligible_results, sep="\t", index=False
    )

    official_results, official_elapsed, official_peak_rss_kb = _run_official(
        eligible,
        rfu_dir=rfu_dir,
        outdir=outdir,
        threshold=threshold,
        rscript_bin=rscript_bin,
    )
    comparison, mismatches, summary = compare(
        scrfu_eligible_results,
        official_results,
        id_col="input_row_id",
        score_tolerance=score_tolerance,
    )
    comparison.to_csv(outdir / "row_comparison.tsv.gz", sep="\t", index=False)
    mismatches.to_csv(outdir / "mismatches.tsv.gz", sep="\t", index=False)

    eligible_scrfu = scrfu_run.df.loc[scrfu_run.df["eligibility_status"].eq("eligible")]
    sequence_counts = eligible_scrfu.groupby("cdr3aa", sort=False)["unique_sequence_id"].nunique()
    reconstruction_mismatches = int(sequence_counts.ne(1).sum())
    order_mismatches = int(
        scrfu_run.df["input_row_id"].astype(str).ne(features["input_row_id"].astype(str)).sum()
    )
    summary.update(
        {
            "mode": "real_execution",
            "total_rows": len(features),
            "unique_sequences": int(eligible["cdr3aa"].nunique()),
            "eligible_sequences": len(eligible),
            "ineligible_sequences": int((~eligible_mask).sum()),
            "below_threshold_sequences": int((~eligible_scrfu["pass_thr"].fillna(False)).sum()),
            "reconstruction_mismatches": reconstruction_mismatches,
            "input_order_reconstruction_mismatches": order_mismatches,
            "scrfu_runtime_seconds": scrfu_elapsed,
            "official_runtime_seconds": official_elapsed,
            "total_runtime_seconds": scrfu_elapsed + official_elapsed,
            "peak_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            "official_child_peak_rss_kb": official_peak_rss_kb,
            "random_seed": seed,
            "threshold": threshold,
            "official_repository_url": _git_value(rfu_dir, "remote", "get-url", "origin"),
            "official_git_commit": _git_value(rfu_dir, "rev-parse", "HEAD"),
            "official_retrieval_date": date.today().isoformat(),
            "official_artifacts": artifacts,
        }
    )
    summary["passed"] = bool(
        summary["passed"] and not reconstruction_mismatches and not order_mismatches
    )
    (outdir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    figure_summary = pd.DataFrame(
        [
            {
                "dataset_label": "adversarial_original_rfu_parity_fixture",
                "input_sha256": _sha256(outdir / "parity_fixture.tsv"),
                "rfu_reference_identifier": artifacts["km5000noMax.Rdata"]["sha256"],
                "software_version": backend.provenance_dict()["scrfu_version"],
                "parameters": json.dumps(
                    {"threshold": threshold, "score_tolerance": score_tolerance}, sort_keys=True
                ),
                "random_seed": seed,
                "total_rows": len(features),
                "eligible_rows": len(eligible),
                "unique_cdr3_queries": int(eligible["cdr3aa"].nunique()),
                "exact_rfu_id_matches": summary["exact_rfu_id_matches"],
                "maximum_score_absolute_error": summary["maximum_score_absolute_error"],
                "threshold_mismatches": summary["threshold_mismatch_count"],
                "reconstruction_mismatches": reconstruction_mismatches + order_mismatches,
                "runtime_seconds": summary["total_runtime_seconds"],
                "peak_rss_kb": summary["peak_rss_kb"],
                "passed": summary["passed"],
            }
        ]
    )
    figure_summary.to_csv(outdir / "figure1_original_rfu_parity.tsv", sep="\t", index=False)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Execute scRFU and untouched original RFU on an adversarial fixture, or compare "
            "two precomputed row-level outputs."
        )
    )
    parser.add_argument("--scrfu-results", type=Path)
    parser.add_argument("--original-results", type=Path)
    parser.add_argument("--rfu-dir", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--id-col", default="input_row_id")
    parser.add_argument("--score-tolerance", type=float, default=1e-12)
    parser.add_argument("--threshold", type=float, default=0.6)
    parser.add_argument("--seed", type=int, default=20260824)
    parser.add_argument("--rscript-bin", default="Rscript")
    args = parser.parse_args()
    if args.score_tolerance < 0:
        parser.error("--score-tolerance must be non-negative")
    if args.rfu_dir is not None:
        if args.scrfu_results is not None or args.original_results is not None:
            parser.error("--rfu-dir cannot be combined with precomputed result paths")
        summary = run_real_parity(
            rfu_dir=args.rfu_dir,
            outdir=args.outdir,
            threshold=args.threshold,
            score_tolerance=args.score_tolerance,
            seed=args.seed,
            rscript_bin=args.rscript_bin,
        )
    else:
        if args.scrfu_results is None or args.original_results is None:
            parser.error("pass --rfu-dir, or both --scrfu-results and --original-results")
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
        summary["mode"] = "precomputed_comparison"
        summary["elapsed_seconds"] = time.perf_counter() - started
        summary["peak_rss_kb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        (args.outdir / "summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n"
        )
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
