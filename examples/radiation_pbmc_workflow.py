from __future__ import annotations

import argparse
import os
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any


def _add_repo_src_to_path() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    src_dir = repo_root / "src"
    if src_dir.exists() and str(src_dir) not in sys.path:
        sys.path.insert(0, str(src_dir))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Template scRFU workflow for user-provided radiation-associated PBMC "
            "single-cell RNA/TCR datasets. This script does not bundle data or "
            "claim biological findings."
        )
    )
    parser.add_argument("--input", required=True, type=Path, help="Input AnnData .h5ad file.")
    parser.add_argument(
        "--rfu-dir",
        required=True,
        type=Path,
        help="Path to an external upstream RFU repository checkout.",
    )
    parser.add_argument(
        "--sample-key",
        default="sample",
        help="Column in adata.obs containing sample identifiers.",
    )
    parser.add_argument(
        "--condition-key",
        default="condition",
        help="Column in adata.obs containing condition, time, or dose labels.",
    )
    parser.add_argument(
        "--outdir",
        default=Path(".scrfu/radiation_pbmc"),
        type=Path,
        help="Output directory for annotated AnnData, TSVs, RFU work files, and plots.",
    )
    return parser


def _fail(message: str) -> int:
    print(f"error: {message}", file=sys.stderr)
    return 2


def _resolve_existing_file(path: Path, label: str) -> Path:
    resolved = path.expanduser()
    if not resolved.exists():
        raise ValueError(f"{label} file not found: {resolved}")
    if not resolved.is_file():
        raise ValueError(f"{label} path is not a file: {resolved}")
    return resolved.resolve()


def _resolve_existing_dir(path: Path, label: str) -> Path:
    resolved = path.expanduser()
    if not resolved.exists():
        raise ValueError(f"{label} directory not found: {resolved}")
    if not resolved.is_dir():
        raise ValueError(f"{label} path is not a directory: {resolved}")
    return resolved.resolve()


def _prepare_runtime(outdir: Path) -> tuple[Any, Any, Any]:
    mpl_config_dir = outdir / ".mplconfig"
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))

    _add_repo_src_to_path()

    import anndata as ad
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    import scrfu

    return ad, plt, scrfu


def _validate_adata(adata: Any, sample_key: str, condition_key: str) -> None:
    if not hasattr(adata, "obsm") or "airr" not in adata.obsm:
        raise ValueError('input AnnData is missing adata.obsm["airr"]')
    if not hasattr(adata, "obs"):
        raise ValueError("input AnnData is missing adata.obs")
    if sample_key not in adata.obs.columns:
        raise ValueError(f'sample-key column "{sample_key}" is missing from adata.obs')
    if condition_key not in adata.obs.columns:
        raise ValueError(f'condition-key column "{condition_key}" is missing from adata.obs')


def _save_plot(path: Path, plotter: Callable[[], Any], plt: Any) -> None:
    ax = plotter()
    ax.figure.tight_layout()
    ax.figure.savefig(path, dpi=150)
    plt.close(ax.figure)


def _write_group_outputs(adata: Any, key: str, label: str, outdir: Path, scrfu: Any) -> None:
    summary = scrfu.tl.rfu_summary(adata, groupby=key)
    summary.to_csv(outdir / f"rfu_summary_by_{label}.tsv", sep="\t", index=False)
    scrfu.tl.aggregate_rfu(adata, groupby=key)
    scrfu.io.export_rfu_matrix(
        adata,
        groupby=key,
        path=outdir / f"rfu_matrix_by_{label}.tsv",
    )


def run_workflow(args: argparse.Namespace) -> int:
    try:
        input_path = _resolve_existing_file(args.input, "input")
        rfu_dir = _resolve_existing_dir(args.rfu_dir, "RFU_DIR")
    except ValueError as exc:
        return _fail(str(exc))

    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    ad, plt, scrfu = _prepare_runtime(outdir)

    print(f"Reading AnnData: {input_path}")
    adata = ad.read_h5ad(input_path)

    try:
        _validate_adata(adata, args.sample_key, args.condition_key)
    except ValueError as exc:
        return _fail(str(exc))

    print("Calling upstream RFU backend")
    scrfu.tl.call_rfu(
        adata,
        backend="rfu_repo",
        rfu_dir=rfu_dir,
        workdir=outdir / "rfu_work",
    )

    print("Writing sample and condition summaries")
    _write_group_outputs(adata, args.sample_key, "sample", outdir, scrfu)
    _write_group_outputs(adata, args.condition_key, "condition", outdir, scrfu)

    print("Writing plots")
    _save_plot(
        outdir / "rfu_bar_by_sample.png",
        lambda: scrfu.pl.rfu_bar(adata, groupby=args.sample_key),
        plt,
    )
    _save_plot(
        outdir / "rfu_heatmap_by_sample.png",
        lambda: scrfu.pl.rfu_heatmap(adata, groupby=args.sample_key),
        plt,
    )
    _save_plot(
        outdir / "rfu_bar_by_condition.png",
        lambda: scrfu.pl.rfu_bar(adata, groupby=args.condition_key),
        plt,
    )
    _save_plot(
        outdir / "rfu_heatmap_by_condition.png",
        lambda: scrfu.pl.rfu_heatmap(adata, groupby=args.condition_key),
        plt,
    )
    _save_plot(outdir / "rfu_score_hist.png", lambda: scrfu.pl.rfu_score_hist(adata), plt)

    annotated_path = outdir / "annotated.h5ad"
    print(f"Writing annotated AnnData: {annotated_path}")
    adata.write_h5ad(annotated_path)

    print(f"Done. Outputs written to: {outdir}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run_workflow(args)


if __name__ == "__main__":
    raise SystemExit(main())
