#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")

from scrfu import __version__, pl, tl  # noqa: E402


def _read_table(path: Path) -> pd.DataFrame:
    separator = "\t" if ".tsv" in path.name.lower() else ","
    return pd.read_csv(path, sep=separator)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_value(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, np.ndarray):
        return [_json_value(item) for item in value.tolist()]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, (pd.Timestamp, datetime)):
        return value.isoformat()
    if value is pd.NA:
        return None
    return value


def _write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(_json_value(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _save_figure(ax: Any, path: Path) -> None:
    ax.figure.tight_layout()
    ax.figure.savefig(path, dpi=200, bbox_inches="tight")
    import matplotlib.pyplot as plt

    plt.close(ax.figure)


def _join_metadata(
    rows: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    cell_key: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if cell_key not in rows or cell_key not in metadata:
        raise ValueError(f"RFU rows and metadata must both contain cell key {cell_key!r}.")
    if metadata[cell_key].isna().any() or metadata[cell_key].duplicated().any():
        raise ValueError("Metadata cell identifiers must be unique and non-missing.")
    overlap = [column for column in metadata if column != cell_key and column in rows]
    if overlap:
        raise ValueError(f"Metadata columns conflict with RFU row columns: {overlap}")
    result_cells = set(rows[cell_key].dropna().astype(str))
    metadata_cells = set(metadata[cell_key].dropna().astype(str))
    joined = rows.merge(metadata, on=cell_key, how="left", sort=False, validate="many_to_one")
    if len(joined) != len(rows):
        raise RuntimeError("Metadata join changed the RFU row count.")
    return joined, {
        "input_row_count": len(rows),
        "output_row_count": len(joined),
        "unmatched_result_cells": len(result_cells - metadata_cells),
        "unmatched_metadata_cells": len(metadata_cells - result_cells),
    }


def run_analysis(
    *,
    rfu_sequences: Path,
    vdjdb: Path,
    vdjdb_release: str,
    outdir: Path,
    rfu_rows: Path | None = None,
    expected_sha256: str | None = None,
    match_mode: str = "cdr3_v",
    v_gene_mode: str = "strip_allele",
    assignment_policy: str = "nearest",
    ambiguity_policy: str = "fractional",
    n_permutations: int = 1000,
    random_state: int = 0,
    metadata: Path | None = None,
    cell_key: str = "cell_id",
    sample_key: str | None = None,
    donor_key: str | None = None,
    phenotype_key: str | None = None,
    condition_key: str | None = None,
    save_permutation_values: bool = False,
) -> dict[str, Any]:
    """Run a local, version-pinned VDJdb evidence analysis without network access."""
    sequences = _read_table(rfu_sequences)
    rows = _read_table(rfu_rows) if rfu_rows is not None else None
    metadata_frame = _read_table(metadata) if metadata is not None else None
    join_qc: dict[str, Any] | None = None
    if rows is not None and metadata_frame is not None:
        rows, join_qc = _join_metadata(rows, metadata_frame, cell_key=cell_key)

    reference = tl.load_vdjdb_reference(
        vdjdb,
        release_label=vdjdb_release,
        expected_sha256=expected_sha256,
    )
    evidence = tl.annotate_vdjdb(
        sequences,
        reference,
        match_mode=match_mode,
        v_gene_mode=v_gene_mode,
    )
    sequence_summaries = tl.summarize_vdjdb_evidence(sequences, evidence)
    row_summary = None
    if rows is not None:
        row_summary = tl.summarize_vdjdb_evidence(rows, evidence).row_summary
    analysis_input = rows if rows is not None else sequences
    coherence = tl.rfu_antigen_coherence(
        analysis_input,
        evidence,
        assignment_policy=assignment_policy,
        ambiguity_policy=ambiguity_policy,
        sample_key=sample_key if rows is not None else None,
        donor_key=donor_key if rows is not None else None,
    )
    abundance = tl.rfu_antigen_abundance(
        analysis_input,
        evidence,
        assignment_policy=assignment_policy,
        ambiguity_policy=ambiguity_policy,
    )
    global_coherence = tl.global_antigen_coherence(
        sequences,
        evidence,
        assignment_policy=assignment_policy,
        ambiguity_policy=ambiguity_policy,
    )
    comparison = tl.compare_antigen_groupings(
        sequences,
        evidence,
        groupings=("rfu", "trbv", "cdr3_length", "trbv_cdr3_length", "size_matched_random"),
        assignment_policy=assignment_policy,
        ambiguity_policy=ambiguity_policy,
        random_state=random_state,
    )
    context = None
    if rows is not None and any((sample_key, donor_key, phenotype_key, condition_key)):
        context = tl.summarize_antigen_context(
            rows,
            evidence,
            cell_key=cell_key,
            sample_key=sample_key,
            donor_key=donor_key,
            phenotype_key=phenotype_key,
            condition_key=condition_key,
        )

    outdir.mkdir(parents=True, exist_ok=True)
    qc = {"provenance": reference.provenance, "validation": reference.validation}
    _write_json(outdir / "vdjdb_reference_qc.json", qc)
    evidence.to_csv(outdir / "vdjdb_matches_long.tsv.gz", sep="\t", index=False)
    sequence_summaries.sequence_summary.to_csv(
        outdir / "sequence_antigen_summary.tsv.gz", sep="\t", index=False
    )
    if row_summary is not None:
        row_summary.to_csv(outdir / "row_antigen_summary.tsv.gz", sep="\t", index=False)
    coherence.to_csv(outdir / "rfu_antigen_coherence.tsv.gz", sep="\t", index=False)
    abundance.to_csv(outdir / "rfu_antigen_abundance.tsv.gz", sep="\t", index=False)
    comparison.to_csv(outdir / "antigen_grouping_comparison.tsv.gz", sep="\t", index=False)
    _write_json(outdir / "global_antigen_coherence.json", global_coherence)

    if context is not None:
        context.sequence_prevalence.to_csv(
            outdir / "sequence_sample_prevalence.tsv.gz", sep="\t", index=False
        )
        context.rfu_antigen_recurrence.to_csv(
            outdir / "rfu_antigen_recurrence.tsv.gz", sep="\t", index=False
        )
        context.phenotype_abundance.to_csv(
            outdir / "rfu_antigen_phenotype_abundance.tsv.gz", sep="\t", index=False
        )
        context.phenotype_evidence_tiers.to_csv(
            outdir / "antigen_evidence_tiers_by_phenotype.tsv.gz", sep="\t", index=False
        )
        context.phenotype_ambiguity.to_csv(
            outdir / "antigen_ambiguity_by_phenotype.tsv.gz", sep="\t", index=False
        )

    permutation = None
    permutation_summary: dict[str, Any]
    try:
        permutation = tl.rfu_antigen_permutation_test(
            sequences,
            evidence,
            n_permutations=n_permutations,
            random_state=random_state,
            assignment_policy=assignment_policy,
            ambiguity_policy=ambiguity_policy,
        )
        permutation_summary = {
            "status": "ok",
            "observed": permutation.observed,
            "null_mean": permutation.null_mean,
            "null_std": permutation.null_std,
            "empirical_upper_tail_probability": permutation.empirical_upper_tail_probability,
            "z_score": permutation.z_score,
            "parameters": permutation.parameters,
        }
        if save_permutation_values:
            pd.DataFrame({"permutation_value": permutation.permutation_values}).to_csv(
                outdir / "permutation_values.tsv.gz", sep="\t", index=False
            )
    except ValueError as exc:
        permutation_summary = {"status": "skipped", "reason": str(exc)}
    _write_json(outdir / "permutation_summary.json", permutation_summary)

    plots: dict[str, str] = {}
    if not abundance.empty:
        _save_figure(pl.rfu_antigen_heatmap(abundance), outdir / "rfu_antigen_heatmap.png")
        _save_figure(pl.rfu_antigen_bubble(abundance), outdir / "rfu_antigen_bubble.png")
        plots["rfu_antigen_heatmap.png"] = "written"
        plots["rfu_antigen_bubble.png"] = "written"
    eligible_coherence = coherence.loc[coherence["eligible_for_coherence"]]
    if not eligible_coherence.empty:
        _save_figure(
            pl.rfu_antigen_coherence(eligible_coherence),
            outdir / "rfu_antigen_coherence.png",
        )
        plots["rfu_antigen_coherence.png"] = "written"
    if permutation is not None:
        _save_figure(
            pl.antigen_permutation_distribution(permutation),
            outdir / "antigen_permutation_distribution.png",
        )
        plots["antigen_permutation_distribution.png"] = "written"

    outputs = sorted(path for path in outdir.iterdir() if path.is_file())
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "scrfu_version": __version__,
        "parameters": {
            "vdjdb_release": vdjdb_release,
            "expected_sha256": expected_sha256,
            "match_mode": match_mode,
            "v_gene_mode": v_gene_mode,
            "assignment_policy": assignment_policy,
            "ambiguity_policy": ambiguity_policy,
            "n_permutations": n_permutations,
            "random_state": random_state,
            "cell_key": cell_key,
            "sample_key": sample_key,
            "donor_key": donor_key,
            "phenotype_key": phenotype_key,
            "condition_key": condition_key,
        },
        "reference": reference.provenance,
        "inputs": {
            "rfu_sequences": str(rfu_sequences.resolve()),
            "rfu_rows": str(rfu_rows.resolve()) if rfu_rows is not None else None,
            "vdjdb": str(vdjdb.resolve()),
            "metadata": str(metadata.resolve()) if metadata is not None else None,
            "rfu_sequence_rows": len(sequences),
            "rfu_result_rows": len(rows) if rows is not None else None,
        },
        "join_qc": join_qc,
        "dimensions": {
            "evidence_rows": len(evidence),
            "matched_unique_sequences": int(evidence["unique_sequence_id"].nunique()),
            "matched_rfus": int(
                sequences.loc[
                    sequences["unique_sequence_id"].isin(evidence["unique_sequence_id"]),
                    "rfu_label",
                ].nunique()
            ),
            "coherence_rows": len(coherence),
            "comparison_rows": len(comparison),
        },
        "permutation": permutation_summary,
        "plots": plots,
        "outputs": {
            path.name: {"sha256": _sha256(path), "size_bytes": path.stat().st_size}
            for path in outputs
        },
    }
    _write_json(outdir / "run_manifest.json", manifest)
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Offline, version-pinned VDJdb evidence analysis for scRFU results."
    )
    parser.add_argument("--rfu-sequences", type=Path, required=True)
    parser.add_argument("--rfu-rows", type=Path)
    parser.add_argument("--vdjdb", type=Path, required=True)
    parser.add_argument("--vdjdb-release", required=True)
    parser.add_argument("--expected-sha256")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--match-mode", choices=("cdr3", "cdr3_v", "paired_exact"), default="cdr3_v"
    )
    parser.add_argument("--v-gene-mode", choices=("exact", "strip_allele"), default="strip_allele")
    parser.add_argument(
        "--assignment-policy", choices=("nearest", "threshold_pass"), default="nearest"
    )
    parser.add_argument(
        "--ambiguity-policy",
        choices=("fractional", "exclude_ambiguous", "multi_label"),
        default="fractional",
    )
    parser.add_argument("--n-permutations", type=int, default=1000)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--save-permutation-values", action="store_true")
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--cell-key", default="cell_id")
    parser.add_argument("--sample-key")
    parser.add_argument("--donor-key")
    parser.add_argument("--phenotype-key")
    parser.add_argument("--condition-key")
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    run_analysis(
        rfu_sequences=args.rfu_sequences,
        rfu_rows=args.rfu_rows,
        vdjdb=args.vdjdb,
        vdjdb_release=args.vdjdb_release,
        expected_sha256=args.expected_sha256,
        outdir=args.outdir,
        match_mode=args.match_mode,
        v_gene_mode=args.v_gene_mode,
        assignment_policy=args.assignment_policy,
        ambiguity_policy=args.ambiguity_policy,
        n_permutations=args.n_permutations,
        random_state=args.random_state,
        save_permutation_values=args.save_permutation_values,
        metadata=args.metadata,
        cell_key=args.cell_key,
        sample_key=args.sample_key,
        donor_key=args.donor_key,
        phenotype_key=args.phenotype_key,
        condition_key=args.condition_key,
    )


if __name__ == "__main__":
    main()
