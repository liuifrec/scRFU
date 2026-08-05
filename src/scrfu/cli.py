from __future__ import annotations

import argparse

from .io import read_h5ad, write_h5ad
from .tl import call_rfu


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

    return p


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
            chain=args.chain,
            airr_key=args.airr_key,
            out_key=args.out_key,
            wrapper_r_path=args.wrapper_r_path,
            rscript_bin=args.rscript_bin,
            extra_r_args=args.extra_r_arg,
            workdir=args.workdir,
        )
        write_h5ad(adata, args.output)
    else:
        raise SystemExit(f"Unknown command: {args.cmd}")
