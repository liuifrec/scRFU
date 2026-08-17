from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


def _add_repo_src_to_path() -> None:
    source = Path(__file__).resolve().parents[1] / "src"
    if source.exists() and str(source) not in sys.path:
        sys.path.insert(0, str(source))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run dataset-independent RFU calling from a receptor cache, table, or H5AD."
    )
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--cache", type=Path, help="Canonical receptor-cache directory")
    inputs.add_argument("--receptors", type=Path, help="Canonical receptor TSV/CSV")
    inputs.add_argument("--input", type=Path, help="H5AD prepared through --adapter")
    parser.add_argument("--adapter", choices=("wells_tcr_ir", "scirpy_airr"))
    parser.add_argument("--airr-key", default="airr")
    parser.add_argument("--tcr-key", default="TCR_IR")
    parser.add_argument("--metadata-column", action="append", default=[])
    parser.add_argument("--chain", default="TRB")
    policy = parser.add_mutually_exclusive_group()
    policy.add_argument("--primary-chain", action="store_true", dest="primary_chain", default=True)
    policy.add_argument("--all-chains", action="store_false", dest="primary_chain")
    parser.add_argument("--rfu-dir", type=Path, default=None)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=20_000)
    parser.add_argument("--threshold", type=float, default=0.6)
    parser.add_argument("--mode", choices=("standard", "map_aware"), default="standard")
    parser.add_argument("--no-resume", action="store_false", dest="resume")
    parser.add_argument("--force-recompute", action="store_true")
    parser.add_argument(
        "--skip-rfu", action="store_true", help="Prepare deterministic unassigned outputs offline"
    )
    parser.add_argument("--groupby", action="append", default=[])
    parser.add_argument("--weighting", choices=("cell", "unique_sequence"), default="cell")
    parser.add_argument(
        "--assignment-policy", choices=("nearest", "threshold_pass"), default="nearest"
    )
    return parser


def _read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t" if ".tsv" in path.name.lower() else ",")


def _offline_results(
    receptors: pd.DataFrame, chain: str
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    selected = receptors.loc[receptors["chain"].astype(str).str.upper().eq(chain.upper())].copy()
    eligible = selected["cdr3aa"].astype("string").str.startswith("C", na=False)
    selected["eligibility_status"] = "ineligible_cdr3_not_starting_c"
    selected.loc[eligible, "eligibility_status"] = "eligible"
    codes, _ = pd.factorize(selected.loc[eligible, "cdr3aa"], sort=False)
    selected["unique_sequence_id"] = pd.Series(pd.NA, index=selected.index, dtype="string")
    lookup = dict(
        zip(
            selected.loc[eligible, "cdr3aa"],
            [f"sequence_{code:08d}" for code in codes],
            strict=True,
        )
    )
    selected.loc[eligible, "unique_sequence_id"] = selected.loc[eligible, "cdr3aa"].map(lookup)
    selected["rfu_id"] = pd.NA
    selected["rfu_label"] = pd.NA
    selected["rfu_score"] = pd.NA
    selected["pass_thr"] = pd.Series(pd.NA, index=selected.index, dtype="boolean")
    selected["rfu_status"] = selected["eligibility_status"].where(~eligible, "skipped_backend")
    mapping = selected.loc[
        :,
        [
            "input_row_id",
            "unique_sequence_id",
            "cell_id",
            "chain",
            "cdr3aa",
            "v_call",
            "source_adapter",
            "source_row_id",
            "eligibility_status",
        ],
    ]
    eligible_rows = selected.loc[eligible]
    per_sequence = eligible_rows.groupby("unique_sequence_id", sort=False, as_index=False).agg(
        cdr3aa=("cdr3aa", "first"),
        query_v_call=("v_call", "first"),
        multiplicity=("input_row_id", "size"),
        rfu_id=("rfu_id", "first"),
        rfu_label=("rfu_label", "first"),
        rfu_score=("rfu_score", "first"),
        pass_thr=("pass_thr", "first"),
        rfu_status=("rfu_status", "first"),
    )
    return per_sequence, selected.reset_index(drop=True), mapping.reset_index(drop=True)


def run(args: argparse.Namespace) -> int:
    _add_repo_src_to_path()
    from scrfu.adapters import prepare_receptors
    from scrfu.io import read_h5ad_dataframe, read_h5ad_obs, read_receptor_cache
    from scrfu.pp import canonicalize_receptor_table, validate_receptor_table
    from scrfu.tl import call_rfu_table, rfu_metrics

    adapter_qc: dict[str, Any] | None = None
    cache_schema_version: int | None = None
    resolved_adapter = args.adapter
    metadata = pd.DataFrame(index=pd.Index([], name="cell_id"))
    if args.cache is not None:
        cached = read_receptor_cache(args.cache)
        receptors, metadata = cached.receptors, cached.cell_metadata
        input_kind = "receptor_cache"
        cache_schema_version = int(cached.manifest["cache_schema_version"])
        resolved_adapter = str(cached.manifest["source_adapter"])
    elif args.receptors is not None:
        receptors = _read_table(args.receptors)
        input_kind = "canonical_table"
    else:
        if args.adapter is None:
            raise ValueError("--adapter is required with --input.")
        if args.adapter == "wells_tcr_ir":
            result = prepare_receptors(
                args.input,
                adapter=args.adapter,
                key=args.tcr_key,
                chain=args.chain,
                primary_chain=args.primary_chain,
                metadata_columns=args.metadata_column,
            )
        else:
            airr = read_h5ad_dataframe(args.input, location="obsm", key=args.airr_key)
            obs = read_h5ad_obs(args.input, columns=args.metadata_column)
            result = prepare_receptors(
                airr,
                adapter=args.adapter,
                airr_key=args.airr_key,
                chain=args.chain,
                primary_chain=args.primary_chain,
                metadata=obs,
                metadata_columns=args.metadata_column,
            )
        receptors, metadata, adapter_qc = result.receptors, result.cell_metadata, result.qc
        input_kind = f"h5ad_{args.adapter}"
    receptors = canonicalize_receptor_table(receptors)
    validate_receptor_table(receptors, strict=True)

    args.outdir.mkdir(parents=True, exist_ok=True)
    receptors.to_csv(args.outdir / "receptors.tsv.gz", sep="\t", index=False, compression="gzip")
    if adapter_qc is not None:
        (args.outdir / "adapter_qc.json").write_text(
            json.dumps(adapter_qc, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    if args.skip_rfu:
        per_sequence, per_row, mapping = _offline_results(receptors, args.chain)
        provenance: dict[str, Any] = {"backend": "skipped", "chain": args.chain}
    else:
        wrapper = Path(__file__).resolve().parents[1] / "r" / "run_rfu_repo.R"
        result = call_rfu_table(
            receptors,
            rfu_dir=args.rfu_dir,
            chain=args.chain,
            mode=args.mode,
            threshold=args.threshold,
            chunk_size=args.chunk_size,
            resume=args.resume,
            force_recompute=args.force_recompute,
            workdir=args.outdir / "backend",
            wrapper_r_path=wrapper,
        )
        per_sequence, per_row, mapping, provenance = (
            result.per_sequence,
            result.per_row,
            result.mapping,
            result.provenance,
        )
    mapping.to_csv(
        args.outdir / "unique_sequence_map.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    per_sequence.to_csv(
        args.outdir / "rfu_results_per_sequence.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    per_row.to_csv(
        args.outdir / "rfu_results_per_row.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    if args.groupby:
        analysis = per_row.merge(
            metadata.reset_index(), on="cell_id", how="left", validate="many_to_one"
        )
        metrics = rfu_metrics(
            analysis,
            groupby=args.groupby,
            weighting=args.weighting,
            chain=args.chain,
            assignment_policy=args.assignment_policy,
        )
        metrics.to_csv(
            args.outdir / "rfu_metrics.tsv.gz", sep="\t", index=False, compression="gzip"
        )
    manifest = {
        **provenance,
        "workflow": "receptor_table",
        "input_kind": input_kind,
        "receptor_row_count": len(receptors),
        "adapter": resolved_adapter,
        "cache_schema_version": cache_schema_version,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "outputs": sorted(path.name for path in args.outdir.iterdir() if path.is_file()),
    }
    (args.outdir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8"
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    try:
        return run(build_parser().parse_args(argv))
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
