from __future__ import annotations

import argparse
from pathlib import Path

from .io import read_h5ad, write_h5ad
from .tl import call_rfu


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="scrfu", description="scRFU: RFU calling for scirpy/AnnData.")
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("call-rfu", help="Run RFU calling and attach results to AnnData.")
    c.add_argument("input", type=str, help="Input .h5ad")
    c.add_argument("-o", "--output", type=str, required=True, help="Output .h5ad")
    c.add_argument("--chain", type=str, default="TRB", help="Chain/locus (default: TRB)")
    c.add_argument("--airr-key", type=str, default="airr", help="obsm key for AIRR table")
    c.add_argument("--out-key", type=str, default="rfu", help="Output key prefix (used for temp files)")
    c.add_argument("--rscript", type=str, required=True, help="Path to RFU_call.R")
    c.add_argument("--atlas", type=str, required=True, help="Path to km5000 centroid atlas TSV.GZ")
    c.add_argument("--extra-r-arg", action="append", default=[], help="Extra arg passed to R script (repeatable)")

    return p


def main(argv: list[str] | None = None) -> None:
    p = build_parser()
    args = p.parse_args(argv)

    if args.cmd == "call-rfu":
        adata = read_h5ad(args.input)
        call_rfu(
            adata,
            chain=args.chain,
            airr_key=args.airr_key,
            out_key=args.out_key,
            rscript_path=Path(args.rscript),
            atlas_path=Path(args.atlas),
            extra_r_args=args.extra_r_arg,
        )
        write_h5ad(adata, args.output)
    else:
        raise SystemExit(f"Unknown command: {args.cmd}")