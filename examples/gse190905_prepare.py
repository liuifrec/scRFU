from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd

CELL_ALIASES = (
    "cell_id",
    "cell",
    "barcode",
    "cellid",
    "cell_barcode",
    "barcodes",
    "cellid",
)
CHAIN_ALIASES = ("chain", "locus")
CDR3_ALIASES = ("cdr3aa", "junction_aa", "cdr3_aa", "cdr3", "cdr3_amino_acid")
V_ALIASES = ("v_call", "v_gene", "trbv", "v", "vregion", "v_region")
PRODUCTIVE_ALIASES = ("productive", "is_productive")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare a user-downloaded GSE190905-like TCR CSV and metadata CSV "
            "for scRFU. No data are downloaded."
        )
    )
    parser.add_argument("--tcr", required=True, type=Path, help="TCR CSV or CSV.GZ file.")
    parser.add_argument("--metadata", required=True, type=Path, help="Metadata CSV or CSV.GZ file.")
    parser.add_argument("--out", required=True, type=Path, help="Output AnnData .h5ad path.")
    parser.add_argument(
        "--report", required=True, type=Path, help="Output preparation report JSON."
    )
    parser.add_argument("--cell-col", default=None, help="TCR cell/barcode column override.")
    parser.add_argument("--chain-col", default=None, help="TCR chain/locus column override.")
    parser.add_argument("--cdr3-col", default=None, help="TCR CDR3 amino-acid column override.")
    parser.add_argument("--v-col", default=None, help="TCR V gene column override.")
    parser.add_argument("--productive-col", default=None, help="TCR productive column override.")
    parser.add_argument(
        "--metadata-cell-col", default=None, help="Metadata cell/barcode column override."
    )
    return parser


def _fail(message: str) -> int:
    print(f"error: {message}", file=sys.stderr)
    return 2


def _normalized(name: object) -> str:
    return "".join(ch for ch in str(name).lower() if ch.isalnum())


def _read_csv(path: Path, label: str) -> pd.DataFrame:
    resolved = path.expanduser()
    if not resolved.exists():
        raise ValueError(f"{label} file not found: {resolved}")
    if not resolved.is_file():
        raise ValueError(f"{label} path is not a file: {resolved}")
    return pd.read_csv(resolved, compression="infer")


def _pick_column(
    df: pd.DataFrame,
    *,
    role: str,
    aliases: tuple[str, ...],
    explicit: str | None,
    source: str,
) -> str:
    if explicit:
        if explicit not in df.columns:
            raise ValueError(f"{source} column override for {role} not found: {explicit}")
        return explicit

    by_normalized = {_normalized(col): col for col in df.columns}
    for alias in aliases:
        col = by_normalized.get(_normalized(alias))
        if col is not None:
            return col
    raise ValueError(
        f"could not infer required {source} column for {role}; "
        f"use an explicit column override. Columns seen: {list(df.columns)}"
    )


def _pick_optional_column(
    df: pd.DataFrame,
    *,
    role: str,
    aliases: tuple[str, ...],
    explicit: str | None,
    source: str,
) -> str | None:
    if explicit:
        if explicit not in df.columns:
            raise ValueError(f"{source} column override for {role} not found: {explicit}")
        return explicit

    by_normalized = {_normalized(col): col for col in df.columns}
    for alias in aliases:
        col = by_normalized.get(_normalized(alias))
        if col is not None:
            return col
    return None


def _truthy_productive(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(False)
    return values.astype(str).str.lower().isin(["true", "t", "1", "yes", "productive"])


def _representative_airr(airr_full: pd.DataFrame, obs_index: pd.Index) -> pd.DataFrame:
    ranked = airr_full.copy()
    ranked["_order"] = range(len(ranked))
    ranked["_is_trb"] = ranked["chain"].astype(str).str.upper().eq("TRB")
    if "productive" in ranked.columns:
        ranked["_productive"] = _truthy_productive(ranked["productive"])
    else:
        ranked["_productive"] = False

    ranked = ranked.sort_values(
        ["cell_id", "_is_trb", "_productive", "_order"],
        ascending=[True, False, False, True],
    )
    representative = ranked.drop_duplicates("cell_id", keep="first").set_index("cell_id")

    aligned = pd.DataFrame(index=obs_index)
    aligned["cell_id"] = obs_index.astype(str)
    reindexed = representative.reindex(obs_index.astype(str))
    for col in ["chain", "cdr3aa", "v_call"]:
        if col in reindexed.columns:
            aligned[col] = reindexed[col].fillna("").astype(str).to_numpy()
    if "productive" in reindexed.columns:
        aligned["productive"] = _truthy_productive(reindexed["productive"]).to_numpy()
    return aligned


def prepare(args: argparse.Namespace) -> int:
    try:
        tcr = _read_csv(args.tcr, "TCR")
        metadata = _read_csv(args.metadata, "metadata")

        tcr_cell_col = _pick_column(
            tcr,
            role="cell",
            aliases=CELL_ALIASES,
            explicit=args.cell_col,
            source="TCR",
        )
        chain_col = _pick_column(
            tcr,
            role="chain",
            aliases=CHAIN_ALIASES,
            explicit=args.chain_col,
            source="TCR",
        )
        cdr3_col = _pick_column(
            tcr,
            role="cdr3",
            aliases=CDR3_ALIASES,
            explicit=args.cdr3_col,
            source="TCR",
        )
        v_col = _pick_column(
            tcr,
            role="v",
            aliases=V_ALIASES,
            explicit=args.v_col,
            source="TCR",
        )
        productive_col = _pick_optional_column(
            tcr,
            role="productive",
            aliases=PRODUCTIVE_ALIASES,
            explicit=args.productive_col,
            source="TCR",
        )
        metadata_cell_col = _pick_column(
            metadata,
            role="cell",
            aliases=CELL_ALIASES,
            explicit=args.metadata_cell_col,
            source="metadata",
        )

        meta = metadata.copy()
        meta[metadata_cell_col] = meta[metadata_cell_col].astype(str)
        meta = meta.drop_duplicates(subset=[metadata_cell_col], keep="first")
        meta = meta.set_index(metadata_cell_col, drop=False)
        meta.index = meta.index.astype(str)

        airr_full = pd.DataFrame(
            {
                "cell_id": tcr[tcr_cell_col].astype(str),
                "chain": tcr[chain_col].astype(str),
                "cdr3aa": tcr[cdr3_col].astype(str),
                "v_call": tcr[v_col].astype(str),
            }
        )
        if productive_col is not None:
            airr_full["productive"] = tcr[productive_col]

        tcr_cells = set(airr_full["cell_id"].dropna().astype(str))
        metadata_cells = set(meta.index.astype(str))
        overlap = tcr_cells.intersection(metadata_cells)
        barcode_overlap_rate = len(overlap) / len(tcr_cells) if tcr_cells else 0.0

        n_productive_trb_rows: int | None = None
        if productive_col is not None:
            trb_rows = airr_full["chain"].astype(str).str.upper().eq("TRB")
            productive_rows = _truthy_productive(airr_full["productive"])
            n_productive_trb_rows = int((trb_rows & productive_rows).sum())

        airr = _representative_airr(airr_full, meta.index)

        import anndata as ad

        adata = ad.AnnData(obs=meta)
        adata.obsm["airr"] = airr

        out_path = args.out.expanduser().resolve()
        report_path = args.report.expanduser().resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.parent.mkdir(parents=True, exist_ok=True)

        adata.write_h5ad(out_path)

        inferred_columns: dict[str, Any] = {
            "tcr_cell_col": tcr_cell_col,
            "chain_col": chain_col,
            "cdr3_col": cdr3_col,
            "v_col": v_col,
            "productive_col": productive_col,
            "metadata_cell_col": metadata_cell_col,
        }
        report = {
            "n_metadata_rows": int(len(metadata)),
            "n_tcr_rows": int(len(tcr)),
            "n_cells_with_tcr": int(len(tcr_cells)),
            "n_overlap_cells": int(len(overlap)),
            "barcode_overlap_rate": barcode_overlap_rate,
            "inferred_columns": inferred_columns,
            "n_productive_trb_rows": n_productive_trb_rows,
            "output_path": str(out_path),
        }
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    except ValueError as exc:
        return _fail(str(exc))

    print(f"Wrote AnnData: {args.out}")
    print(f"Wrote report: {args.report}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return prepare(args)


if __name__ == "__main__":
    raise SystemExit(main())
