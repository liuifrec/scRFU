"""Verify representative real VDJdb outputs on bounded deterministic subsets."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from scrfu import __version__, tl


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


def _verify_one(
    config: dict[str, Any], reference: tl.VDJdbReference, *, subset_size: int, seed: int
) -> dict[str, Any]:
    run_dir = Path(config["run_dir"])
    manifest = json.loads((run_dir / "run_manifest.json").read_text())
    evidence_path = run_dir / "vdjdb_matches_long.tsv.gz"
    recorded_hash = manifest["outputs"][evidence_path.name]["sha256"]
    actual_hash = _sha256(evidence_path)
    stored = _read(evidence_path)
    selected_ids = sorted(stored["unique_sequence_id"].astype(str).unique())[:subset_size]
    selected = set(selected_ids)
    rows = _read(Path(config["rows"]))
    rows = rows.loc[rows["unique_sequence_id"].astype(str).isin(selected)].copy()
    sequences = _read(Path(config["sequences"]))
    sequences = sequences.loc[sequences["unique_sequence_id"].astype(str).isin(selected)].copy()
    reproduced = tl.annotate_vdjdb(
        rows,
        reference,
        match_mode=config["match_mode"],
        v_gene_mode="strip_allele",
        expand_rows=False,
    )
    expected = stored.loc[stored["unique_sequence_id"].astype(str).isin(selected)].copy()
    keys = ["unique_sequence_id", "match_query_id", "reference_row_id"]
    expected_keys = (
        expected[keys].astype(str).sort_values(keys, kind="stable").reset_index(drop=True)
    )
    reproduced_keys = (
        reproduced[keys].astype(str).sort_values(keys, kind="stable").reset_index(drop=True)
    )
    identity_exact = expected_keys.equals(reproduced_keys)
    summary = tl.summarize_vdjdb_evidence(rows, reproduced)
    reconstruction_exact = (
        len(summary.row_summary) == len(rows)
        and summary.row_summary["input_row_id"].astype(str).tolist()
        == rows["input_row_id"].astype(str).tolist()
    )
    query_variant_count = int(reproduced["match_query_id"].nunique())
    sequence_identity_count = int(reproduced["unique_sequence_id"].nunique())
    coherence_a = tl.global_antigen_coherence(
        sequences,
        reproduced,
        assignment_policy=config["assignment_policy"],
        ambiguity_policy="fractional",
    )
    coherence_b = tl.global_antigen_coherence(
        sequences,
        reproduced,
        assignment_policy=config["assignment_policy"],
        ambiguity_policy="fractional",
    )
    coherence_exact = coherence_a.keys() == coherence_b.keys() and all(
        (math.isnan(left) and math.isnan(right))
        if isinstance(left, float)
        and isinstance(right, float)
        and (math.isnan(left) or math.isnan(right))
        else left == right
        for left, right in ((coherence_a[key], coherence_b[key]) for key in coherence_a)
    )
    null_exact: bool | None = None
    null_reason: str | None = None
    try:
        null_a = tl.rfu_antigen_permutation_test(
            sequences,
            reproduced,
            n_permutations=25,
            random_state=seed,
            assignment_policy=config["assignment_policy"],
            ambiguity_policy="fractional",
        )
        null_b = tl.rfu_antigen_permutation_test(
            sequences,
            reproduced,
            n_permutations=25,
            random_state=seed,
            assignment_policy=config["assignment_policy"],
            ambiguity_policy="fractional",
        )
        null_exact = bool((null_a.permutation_values == null_b.permutation_values).all())
    except ValueError as error:
        null_reason = str(error)
    return {
        "label": config["label"],
        "match_mode": config["match_mode"],
        "assignment_policy": config["assignment_policy"],
        "stored_manifest_sha256": _sha256(run_dir / "run_manifest.json"),
        "stored_evidence_sha256": actual_hash,
        "stored_evidence_hash_matches_manifest": actual_hash == recorded_hash,
        "subset_sequence_count": len(selected_ids),
        "subset_input_row_count": len(rows),
        "expected_evidence_rows": len(expected),
        "reproduced_evidence_rows": len(reproduced),
        "evidence_identity_exact": identity_exact,
        "row_reconstruction_exact": reconstruction_exact,
        "rfu_sequence_identity_count": sequence_identity_count,
        "match_query_variant_count": query_variant_count,
        "coherence_repeat_exact": coherence_exact,
        "null_seed_repeat_exact": null_exact,
        "null_skip_reason": null_reason,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vdjdb", type=Path, required=True)
    parser.add_argument("--release", required=True)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--wells-root", type=Path, required=True)
    parser.add_argument("--gse190905-root", type=Path, required=True)
    parser.add_argument("--gse157007-root", type=Path, required=True)
    parser.add_argument("--subset-size", type=int, default=64)
    parser.add_argument("--random-state", type=int, default=20260825)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    reference = tl.load_vdjdb_reference(
        args.vdjdb,
        release_label=args.release,
        expected_sha256=args.expected_sha256,
    )
    configs = [
        {
            "label": "Wells CDR3 nearest fractional",
            "run_dir": args.vdjdb.parent.parent / "analysis/wells_full/cdr3__nearest__fractional",
            "rows": args.wells_root / "rfu_results_per_row.tsv.gz",
            "sequences": args.wells_root / "rfu_results_per_sequence.tsv.gz",
            "match_mode": "cdr3",
            "assignment_policy": "nearest",
        },
        {
            "label": "Wells CDR3+V threshold fractional",
            "run_dir": args.vdjdb.parent.parent
            / "analysis/wells_full/cdr3_v__threshold_pass__fractional",
            "rows": args.wells_root / "rfu_results_per_row.tsv.gz",
            "sequences": args.wells_root / "rfu_results_per_sequence.tsv.gz",
            "match_mode": "cdr3_v",
            "assignment_policy": "threshold_pass",
        },
        {
            "label": "GSE190905 CDR3 nearest fractional",
            "run_dir": args.gse190905_root / "vdjdb/cdr3__nearest__fractional",
            "rows": args.gse190905_root / "rfu/rfu_results_per_row.tsv.gz",
            "sequences": args.gse190905_root / "rfu/rfu_results_per_sequence.tsv.gz",
            "match_mode": "cdr3",
            "assignment_policy": "nearest",
        },
        {
            "label": "GSE157007 CDR3 nearest fractional",
            "run_dir": args.gse157007_root / "vdjdb/cdr3__nearest__fractional",
            "rows": args.gse157007_root / "rfu/rfu_results_per_row.tsv.gz",
            "sequences": args.gse157007_root / "rfu/rfu_results_per_sequence.tsv.gz",
            "match_mode": "cdr3",
            "assignment_policy": "nearest",
        },
    ]
    checks = [
        _verify_one(config, reference, subset_size=args.subset_size, seed=args.random_state)
        for config in configs
    ]
    status = (
        "valid"
        if all(
            check["stored_evidence_hash_matches_manifest"]
            and check["evidence_identity_exact"]
            and check["row_reconstruction_exact"]
            and check["coherence_repeat_exact"]
            and check["null_seed_repeat_exact"] is not False
            for check in checks
        )
        else "invalid"
    )
    report = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "scrfu_version": __version__,
        "reference_release": args.release,
        "reference_sha256": reference.provenance["sha256"],
        "subset_size": args.subset_size,
        "random_state": args.random_state,
        "status": status,
        "checks": checks,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    if status != "valid":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
