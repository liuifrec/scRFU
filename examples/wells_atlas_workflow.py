from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pandas as pd


def _add_repo_src_to_path() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    src_dir = repo_root / "src"
    if src_dir.exists() and str(src_dir) not in sys.path:
        sys.path.insert(0, str(src_dir))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run the Wells TCR_IR adapter and restartable scRFU assignment from a "
            "user-supplied public-atlas H5AD or prepared receptor cache."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Input Wells atlas H5AD or prepared receptor-cache directory.",
    )
    parser.add_argument(
        "--source-h5ad",
        type=Path,
        default=None,
        help="Optional source H5AD used to validate a prepared cache fingerprint.",
    )
    parser.add_argument(
        "--obs-column",
        action="append",
        default=[],
        help="Observation metadata column to retain from H5AD input (repeatable).",
    )
    parser.add_argument(
        "--rfu-dir",
        type=Path,
        default=None,
        help="Official-compatible RFU checkout; falls back to RFU_DIR.",
    )
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory.")
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=20_000,
        help="Unique CDR3 queries per serial restartable chunk (default: 20000).",
    )
    chain_group = parser.add_mutually_exclusive_group()
    chain_group.add_argument(
        "--primary-chain",
        action="store_true",
        dest="primary_chain",
        default=True,
        help="Use the primary VDJ_1 chain per cell (default).",
    )
    chain_group.add_argument(
        "--all-productive-chains",
        action="store_false",
        dest="primary_chain",
        help="Extract productive TRB chains from VDJ_1 and VDJ_2 fields.",
    )
    parser.add_argument("--threshold", type=float, default=0.6, help="RFU threshold.")
    parser.add_argument(
        "--mode",
        choices=("standard", "map_aware"),
        default="standard",
        help="Backend mode; standard is canonical and default.",
    )
    parser.add_argument("--max-cells", type=int, default=None, help="Optional smoke-test limit.")
    parser.add_argument(
        "--no-resume",
        action="store_false",
        dest="resume",
        help="Do not reuse validated completed chunks.",
    )
    parser.add_argument(
        "--force-recompute",
        action="store_true",
        help="Recompute all chunks; takes precedence over resume.",
    )
    parser.add_argument(
        "--write-annotated",
        action="store_true",
        help="Write annotated.h5ad (primary-chain mode only).",
    )
    return parser


def _atomic_write_json(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def _subset_for_smoke(adata: Any, max_cells: int | None) -> Any:
    if max_cells is None:
        return adata
    if max_cells <= 0:
        raise ValueError("--max-cells must be a positive integer.")
    selected = list(map(str, adata.obs_names[:max_cells]))
    uns_tcr = None
    if "TCR_IR" in adata.uns:
        value = adata.uns["TCR_IR"]
        if isinstance(value, pd.DataFrame):
            uns_tcr = value.loc[selected].copy()
        elif isinstance(value, Mapping):
            uns_tcr = pd.DataFrame(dict(value)).iloc[: len(selected)].copy()
    subset = adata[: len(selected)].copy()
    if uns_tcr is not None:
        subset.uns["TCR_IR"] = uns_tcr
    return subset


def _subset_receptor_data(data: Any, max_cells: int | None) -> Any:
    if max_cells is None or max_cells >= len(data.obs_names):
        return data
    if max_cells <= 0:
        raise ValueError("--max-cells must be a positive integer.")
    from scrfu.wells import WellsReceptorData

    return WellsReceptorData(
        obs=data.obs.iloc[:max_cells].copy(),
        tcr_ir=data.tcr_ir.iloc[:max_cells].copy(),
        atlas_shape=data.atlas_shape,
        tcr_ir_container=data.tcr_ir_container,
        source_h5ad=data.source_h5ad,
        cache_manifest=data.cache_manifest,
    )


def _tcr_ir_frame(adata: Any) -> tuple[pd.DataFrame, str]:
    if "TCR_IR" in adata.uns:
        value = adata.uns["TCR_IR"]
        source = "uns"
    elif "TCR_IR" in adata.obsm:
        value = adata.obsm["TCR_IR"]
        source = "obsm"
    else:
        raise ValueError('Input AnnData has neither uns["TCR_IR"] nor obsm["TCR_IR"].')
    frame = value.copy() if isinstance(value, pd.DataFrame) else pd.DataFrame(value)
    if len(frame) != len(adata.obs_names):
        raise ValueError(
            f"TCR_IR row count ({len(frame)}) does not match AnnData observations "
            f"({len(adata.obs_names)})."
        )
    if isinstance(frame.index, pd.RangeIndex):
        frame.index = pd.Index(map(str, adata.obs_names))
    elif list(map(str, frame.index)) != list(map(str, adata.obs_names)):
        raise ValueError("TCR_IR index does not exactly match adata.obs_names.")
    return frame, source


def _column(frame: pd.DataFrame, slot: int, field: str) -> Any | None:
    normalized = {str(column).lower(): column for column in frame.columns}
    aliases = [f"tcr-ir_vdj_{slot}_{field}", f"ir_vdj_{slot}_{field}"]
    if field == "junction_aa":
        aliases.extend([f"tcr-ir_vdj_{slot}_cdr3", f"ir_vdj_{slot}_cdr3"])
    for alias in aliases:
        if alias in normalized:
            return normalized[alias]
    return None


def _truthy(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(False)
    return values.astype(str).str.lower().isin(["true", "t", "1", "yes"])


def extract_wells_features(
    adata: Any, *, primary_chain: bool
) -> tuple[pd.DataFrame, dict[str, Any]]:
    from scrfu.adapters import adapt_wells_tcr_ir

    result = adapt_wells_tcr_ir(adata, chain="TRB", primary_chain=primary_chain)
    features = result.receptors.rename(columns={"v_call": "trbv", "source_slot": "chain_slot"}).loc[
        :, ["cell_id", "cdr3aa", "trbv", "chain_slot"]
    ]
    qc = {
        **result.qc,
        "n_obs": len(adata.obs_names),
        "tcr_ir_rows": result.qc["source_row_count"],
        "extracted_trb_rows": len(features),
        "extracted_cells": int(features["cell_id"].nunique()) if not features.empty else 0,
    }
    return features, qc


def run_workflow(args: argparse.Namespace) -> int:
    _add_repo_src_to_path()

    from scrfu.adapters import adapt_wells_tcr_ir
    from scrfu.tl import call_rfu_table
    from scrfu.wells import load_wells_receptor_cache, read_wells_receptors_h5ad

    input_path = args.input.expanduser().resolve()
    if input_path.is_dir():
        source_h5ad = (
            args.source_h5ad.expanduser().resolve() if args.source_h5ad is not None else None
        )
        receptor_data = load_wells_receptor_cache(input_path, source_h5ad=source_h5ad)
        receptor_data = _subset_receptor_data(receptor_data, args.max_cells)
        input_kind = "prepared_wells_cache"
    elif input_path.is_file():
        if args.source_h5ad is not None:
            raise ValueError("--source-h5ad is only valid when --input is a prepared cache.")
        receptor_data = read_wells_receptors_h5ad(
            input_path,
            obs_columns=args.obs_column,
            max_cells=args.max_cells,
        )
        input_kind = "h5ad_targeted_reader"
    else:
        raise ValueError(f"Input H5AD or receptor cache not found: {input_path}")
    if args.write_annotated and input_kind != "h5ad_targeted_reader":
        raise ValueError("--write-annotated requires original H5AD input, not a receptor cache.")
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    adapted = adapt_wells_tcr_ir(receptor_data, chain="TRB", primary_chain=args.primary_chain)
    receptors = adapted.receptors
    adapter_qc = adapted.qc
    features = receptors.rename(columns={"v_call": "trbv", "source_slot": "chain_slot"}).loc[
        :, ["cell_id", "cdr3aa", "trbv", "chain_slot"]
    ]
    if features.empty:
        raise ValueError("Wells TCR_IR adapter extracted no productive TRB rows.")
    if args.write_annotated and not args.primary_chain:
        raise ValueError("--write-annotated currently requires --primary-chain.")

    extracted_path = outdir / "extracted_trb.tsv.gz"
    features.to_csv(extracted_path, sep="\t", index=False, compression="gzip")
    receptors.to_csv(outdir / "receptors.tsv.gz", sep="\t", index=False, compression="gzip")
    _atomic_write_json(outdir / "adapter_qc.json", adapter_qc)

    wrapper = Path(__file__).resolve().parents[1] / "r" / "run_rfu_repo.R"
    table_run = call_rfu_table(
        receptors,
        rfu_dir=args.rfu_dir,
        chain="TRB",
        mode=args.mode,
        threshold=args.threshold,
        deduplicate=True,
        chunk_size=args.chunk_size,
        resume=args.resume,
        force_recompute=args.force_recompute,
        workdir=outdir / "backend",
        wrapper_r_path=wrapper,
    )
    per_cell = table_run.per_row
    eligible = per_cell[per_cell["eligibility_status"].eq("eligible")]
    sequence_map = per_cell.loc[
        :,
        [
            "input_row_id",
            "unique_sequence_id",
            "cell_id",
            "cdr3aa",
            "trbv",
            "eligibility_status",
        ],
    ]
    sequence_summary = eligible.groupby("unique_sequence_id", sort=False, as_index=False).agg(
        cdr3aa=("cdr3aa", "first"), query_trbv=("trbv", "first"), multiplicity=("cell_id", "size")
    )
    sequence_assignments = eligible.drop_duplicates("unique_sequence_id", keep="first").loc[
        :,
        [
            "unique_sequence_id",
            "rfu_id",
            "rfu_label",
            "rfu_score",
            "pass_thr",
            "rfu_status",
        ],
    ]
    sequence_results = sequence_summary.merge(
        sequence_assignments,
        on="unique_sequence_id",
        how="left",
        sort=False,
        validate="one_to_one",
    )
    sequence_map.to_csv(
        outdir / "unique_sequence_map.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    sequence_results.to_csv(
        outdir / "rfu_results_per_sequence.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    per_cell.to_csv(
        outdir / "rfu_results_per_cell.tsv.gz", sep="\t", index=False, compression="gzip"
    )
    per_cell.to_csv(
        outdir / "rfu_results_per_row.tsv.gz", sep="\t", index=False, compression="gzip"
    )

    provenance = {
        **table_run.provenance,
        "workflow": "wells_atlas",
        "input": str(input_path),
        "input_kind": input_kind,
        "source_atlas_dimensions": list(receptor_data.atlas_shape),
        "selected_metadata_columns": list(receptor_data.obs.columns),
        "primary_chain": args.primary_chain,
        "max_cells": args.max_cells,
        "outputs": [
            "extracted_trb.tsv.gz",
            "receptors.tsv.gz",
            "adapter_qc.json",
            "unique_sequence_map.tsv.gz",
            "rfu_results_per_sequence.tsv.gz",
            "rfu_results_per_cell.tsv.gz",
            "rfu_results_per_row.tsv.gz",
        ],
    }
    _atomic_write_json(outdir / "run_manifest.json", provenance)
    if args.write_annotated:
        import anndata as ad

        from scrfu.attach import attach_rfu_results

        # Explicit opt-in: writing annotated output necessarily reads expression data.
        adata = _subset_for_smoke(ad.read_h5ad(input_path), args.max_cells)
        attach_rfu_results(adata, features, per_cell, provenance=provenance)
        if hasattr(ad.settings, "allow_write_nullable_strings"):
            ad.settings.allow_write_nullable_strings = True
        adata.write_h5ad(outdir / "annotated.h5ad")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return run_workflow(args)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
