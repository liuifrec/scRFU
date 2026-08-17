from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .adapters import list_receptor_adapters, prepare_receptors
from .io import (
    migrate_wells_receptor_cache,
    read_h5ad,
    read_h5ad_dataframe,
    read_h5ad_obs,
    write_h5ad,
    write_receptor_cache,
)
from .pp import validate_receptor_table
from .tl import call_rfu
from .wells import prepare_wells_receptor_cache


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="scrfu", description="scRFU: RFU calling for scirpy/AnnData.")
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("call-rfu", help="Run RFU calling and attach results to AnnData.")
    c.add_argument("input", type=str, help="Input .h5ad")
    c.add_argument("-o", "--output", type=str, required=True, help="Output .h5ad")
    c.add_argument(
        "--rfu-dir",
        type=str,
        default=None,
        help="Path to upstream RFU checkout (falls back to RFU_DIR)",
    )
    c.add_argument(
        "--mode",
        choices=("standard", "map_aware"),
        default="standard",
        help="RFU backend mode (default: standard)",
    )
    c.add_argument("--threshold", type=float, default=0.6, help="RFU threshold (default: 0.6)")
    c.add_argument(
        "--no-deduplicate",
        action="store_false",
        dest="deduplicate",
        help="Submit every eligible row instead of unique CDR3 queries",
    )
    c.add_argument(
        "--chunk-size",
        type=int,
        default=None,
        help="Unique eligible CDR3 queries per chunk",
    )
    c.add_argument("--max-workers", type=int, default=1, help="Independent chunk workers")
    c.add_argument(
        "--executor", choices=("process", "thread"), default="process", help="Chunk executor"
    )
    c.add_argument(
        "--no-resume",
        action="store_false",
        dest="resume",
        help="Do not reuse previously completed chunks",
    )
    c.add_argument(
        "--force-recompute",
        action="store_true",
        help="Recompute every chunk; takes precedence over resume",
    )
    c.add_argument("--chain", type=str, default="TRB", help="Chain/locus (default: TRB)")
    c.add_argument("--airr-key", type=str, default="airr", help="obsm key for AIRR table")
    c.add_argument("--out-key", type=str, default="rfu", help="Output provenance key")
    c.add_argument(
        "--wrapper-r-path",
        type=str,
        default="r/run_rfu_repo.R",
        help="Path to scrfu R wrapper script",
    )
    c.add_argument("--rscript-bin", type=str, default="Rscript", help="Rscript executable")
    c.add_argument("--workdir", type=str, default=None, help="Scratch directory for RFU run")
    c.add_argument(
        "--extra-r-arg",
        action="append",
        default=[],
        help="Extra arg passed to R script (repeatable)",
    )

    w = sub.add_parser(
        "prepare-wells",
        help="Extract a compact Wells TCR_IR/metadata cache without reading expression matrices.",
    )
    w.add_argument("input", type=str, help="Source Wells atlas .h5ad")
    w.add_argument("-o", "--output-dir", type=str, required=True, help="Cache output directory")
    w.add_argument(
        "--obs-column",
        action="append",
        default=[],
        help="Observation metadata column to retain (repeatable)",
    )

    generic = sub.add_parser(
        "prepare-receptors",
        help="Prepare a dataset-independent canonical receptor cache with a named adapter.",
    )
    generic.add_argument("--input", type=str, default=None, help="Input H5AD or AIRR table")
    generic.add_argument("--adapter", type=str, default=None, help="Named receptor adapter")
    generic.add_argument("-o", "--outdir", type=str, default=None, help="Cache output directory")
    generic.add_argument("--chain", type=str, default="TRB", help="Chain/locus selection")
    policy = generic.add_mutually_exclusive_group()
    policy.add_argument("--primary-chain", action="store_true", dest="primary_chain", default=True)
    policy.add_argument("--all-chains", action="store_false", dest="primary_chain")
    generic.add_argument(
        "--metadata-columns", nargs="*", default=[], help="Observation metadata columns to retain"
    )
    generic.add_argument("--airr-key", default="airr", help="obsm key for AIRR data")
    generic.add_argument("--tcr-key", default="TCR_IR", help="Wells TCR_IR key")
    generic.add_argument(
        "--include-non-cells",
        action="store_false",
        dest="filter_is_cell",
        default=True,
        help="For Cell Ranger, retain rows where is_cell is false",
    )
    generic.add_argument(
        "--include-low-confidence",
        action="store_false",
        dest="filter_high_confidence",
        default=True,
        help="For Cell Ranger, retain rows where high_confidence is false",
    )
    generic.add_argument(
        "--include-nonproductive",
        action="store_false",
        dest="productive_only",
        default=True,
        help="Retain nonproductive receptor rows",
    )
    generic.add_argument("--force", action="store_true", help="Replace an existing cache")
    generic.add_argument(
        "--list-adapters", action="store_true", help="List built-in adapters and exit"
    )
    generic.add_argument(
        "--validate-only", action="store_true", help="Validate preparation without writing a cache"
    )
    generic.add_argument("--max-cells", type=int, default=None, help="Optional H5AD prefix limit")

    migrate = sub.add_parser(
        "migrate-receptor-cache", help="Explicitly migrate a legacy Wells cache."
    )
    migrate.add_argument("input", help="Legacy Wells cache directory")
    migrate.add_argument("-o", "--outdir", required=True, help="New receptor cache directory")
    migrate.add_argument("--force", action="store_true")
    w.add_argument(
        "--max-cells",
        type=int,
        default=None,
        help="Optional prefix limit for a smoke-test cache",
    )

    return p


def _prepare_receptors_command(args: argparse.Namespace) -> dict[str, object]:
    if args.list_adapters:
        for name in list_receptor_adapters():
            print(name)
        return {"adapters": list_receptor_adapters()}
    if args.input is None or args.adapter is None:
        raise ValueError("--input and --adapter are required unless --list-adapters is used.")
    if args.outdir is None and not args.validate_only:
        raise ValueError("--outdir is required unless --validate-only is used.")
    path = Path(args.input).expanduser().resolve()
    adapter = args.adapter.strip().lower()
    source: object
    kwargs: dict[str, object] = {
        "chain": args.chain,
        "primary_chain": args.primary_chain,
        "metadata_columns": args.metadata_columns,
    }
    if adapter in {"wells", "wells_tcr_ir"}:
        source = path
        kwargs["key"] = args.tcr_key
        if args.max_cells is not None:
            from .wells import read_wells_receptors_h5ad

            source = read_wells_receptors_h5ad(
                path, obs_columns=args.metadata_columns, max_cells=args.max_cells
            )
    elif adapter in {"cellranger_vdj", "cellranger", "tenx_vdj"}:
        source = path
        kwargs.update(
            {
                "filter_is_cell": args.filter_is_cell,
                "filter_high_confidence": args.filter_high_confidence,
                "productive_only": args.productive_only,
            }
        )
    elif path.suffix.lower() == ".h5ad":
        airr = read_h5ad_dataframe(
            path, location="obsm", key=args.airr_key, max_rows=args.max_cells
        )
        obs = read_h5ad_obs(path, columns=args.metadata_columns, max_rows=args.max_cells)
        source = airr
        kwargs.update({"airr_key": args.airr_key, "metadata": obs})
    else:
        separator = "\t" if ".tsv" in path.name.lower() else ","
        source = pd.read_csv(path, sep=separator)
        kwargs["airr_key"] = args.airr_key
    result = prepare_receptors(source, adapter=adapter, **kwargs)
    report = validate_receptor_table(result.receptors, strict=True)
    summary: dict[str, object] = {
        "adapter": result.adapter_name,
        "receptor_rows": len(result.receptors),
        "unique_cells": int(result.receptors["cell_id"].nunique()),
        "validation": report["status"],
    }
    if not args.validate_only:
        provenance = result.provenance
        manifest = write_receptor_cache(
            args.outdir,
            result.receptors,
            result.cell_metadata,
            source_adapter=result.adapter_name,
            source_adapter_version=result.adapter_version,
            source_format=str(provenance.get("source_format", "unknown")),
            source_path=path,
            adapter_qc=result.qc,
            adapter_configuration={
                "chain": args.chain,
                "primary_chain": args.primary_chain,
                "airr_key": args.airr_key,
                "tcr_key": args.tcr_key,
                "filter_is_cell": args.filter_is_cell,
                "filter_high_confidence": args.filter_high_confidence,
                "productive_only": args.productive_only,
            },
            selected_metadata_columns=args.metadata_columns,
            source_atlas_dimensions=provenance.get("source_atlas_dimensions"),
            force=args.force,
        )
        summary["cache_schema_version"] = manifest["cache_schema_version"]
    print(json.dumps(summary, sort_keys=True))
    return summary


def main(argv: list[str] | None = None) -> None:
    p = build_parser()
    args = p.parse_args(argv)

    if args.cmd == "call-rfu":
        adata = read_h5ad(args.input)
        call_rfu(
            adata,
            rfu_dir=args.rfu_dir,
            mode=args.mode,
            threshold=args.threshold,
            deduplicate=args.deduplicate,
            chunk_size=args.chunk_size,
            max_workers=args.max_workers,
            executor=args.executor,
            resume=args.resume,
            force_recompute=args.force_recompute,
            chain=args.chain,
            airr_key=args.airr_key,
            out_key=args.out_key,
            wrapper_r_path=args.wrapper_r_path,
            rscript_bin=args.rscript_bin,
            extra_r_args=args.extra_r_arg,
            workdir=args.workdir,
        )
        write_h5ad(adata, args.output)
    elif args.cmd == "prepare-wells":
        prepare_wells_receptor_cache(
            args.input,
            args.output_dir,
            obs_columns=args.obs_column,
            max_cells=args.max_cells,
        )
    elif args.cmd == "prepare-receptors":
        _prepare_receptors_command(args)
    elif args.cmd == "migrate-receptor-cache":
        migrate_wells_receptor_cache(args.input, args.outdir, force=args.force)
    else:
        raise SystemExit(f"Unknown command: {args.cmd}")


if __name__ == "__main__":
    main()
