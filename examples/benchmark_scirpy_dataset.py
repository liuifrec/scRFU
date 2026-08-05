from __future__ import annotations

import argparse
import json
import os
import sys
from collections.abc import Callable
from datetime import datetime, timezone
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
            "Run a reproducible scRFU benchmark workflow on a user-provided "
            "AnnData/scirpy-style h5ad file. No data are downloaded."
        )
    )
    parser.add_argument("--input", required=True, type=Path, help="Input AnnData .h5ad file.")
    parser.add_argument(
        "--rfu-dir",
        required=True,
        type=Path,
        help="Path to an external upstream RFU repository checkout.",
    )
    parser.add_argument("--groupby", required=True, help="adata.obs column for grouped outputs.")
    parser.add_argument(
        "--cell-type-key",
        default=None,
        help="Optional adata.obs column for cell-type outputs.",
    )
    parser.add_argument("--airr-key", default="airr", help="adata.obsm key containing AIRR data.")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory.")
    parser.add_argument(
        "--skip-rfu",
        action="store_true",
        help="Skip RFU execution and use existing adata.obs rfu_label/rfu_score columns.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Load input, run AIRR validation, write manifest, and skip RFU/plot outputs.",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=None,
        help="Optional maximum number of cells for quick testing.",
    )
    parser.add_argument(
        "--write-annotated",
        action="store_true",
        help="Write annotated.h5ad after running the benchmark.",
    )
    return parser


def _fail(message: str) -> int:
    print(f"error: {message}", file=sys.stderr)
    return 2


def _resolve_input(path: Path) -> Path:
    resolved = path.expanduser()
    if not resolved.exists():
        raise ValueError(f"input file not found: {resolved}")
    if not resolved.is_file():
        raise ValueError(f"input path is not a file: {resolved}")
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


def _save_plot(path: Path, plotter: Callable[[], Any], plt: Any) -> None:
    ax = plotter()
    ax.figure.tight_layout()
    ax.figure.savefig(path, dpi=150)
    plt.close(ax.figure)


def _require_obs_columns(adata: Any, columns: list[str]) -> None:
    missing = [col for col in columns if col not in adata.obs.columns]
    if missing:
        raise ValueError(f"adata.obs is missing required columns: {missing}")


def _write_manifest(
    *,
    path: Path,
    scrfu: Any,
    input_path: Path,
    outdir: Path,
    rfu_dir: Path,
    rfu_run: bool,
    groupby: str,
    cell_type_key: str | None,
    n_cells: int,
    output_files: list[Path],
    notes: list[str],
) -> None:
    manifest = {
        "scrfu_version": scrfu.__version__,
        "input_path": str(input_path),
        "output_directory": str(outdir),
        "rfu_dir": str(rfu_dir.expanduser()),
        "rfu_run": rfu_run,
        "groupby": groupby,
        "cell_type_key": cell_type_key,
        "n_cells": n_cells,
        "output_files": [str(p) for p in output_files],
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "notes": [
            "scRFU depends on a user-provided upstream RFU checkout; no RFU code or data are bundled.",
            *notes,
        ],
    }
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def _write_group_outputs(
    adata: Any,
    *,
    groupby: str,
    label: str,
    outdir: Path,
    scrfu: Any,
    plt: Any,
) -> list[Path]:
    outputs: list[Path] = []

    summary_path = outdir / f"rfu_summary_by_{label}.tsv"
    matrix_path = outdir / f"rfu_matrix_by_{label}.tsv"
    heatmap_path = outdir / f"rfu_heatmap_by_{label}.png"

    scrfu.tl.rfu_summary(adata, groupby=groupby).to_csv(summary_path, sep="\t", index=False)
    scrfu.io.export_rfu_matrix(adata, groupby=groupby, path=matrix_path)
    _save_plot(heatmap_path, lambda: scrfu.pl.rfu_heatmap(adata, groupby=groupby), plt)

    outputs.extend([summary_path, matrix_path, heatmap_path])
    return outputs


def run_benchmark(args: argparse.Namespace) -> int:
    try:
        input_path = _resolve_input(args.input)
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        return _fail(str(exc))

    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    ad, plt, scrfu = _prepare_runtime(outdir)

    try:
        adata = ad.read_h5ad(input_path)
        if args.max_cells is not None:
            if args.max_cells <= 0:
                raise ValueError("--max-cells must be a positive integer")
            adata = adata[: args.max_cells].copy()

        outputs: list[Path] = []
        notes: list[str] = []

        validate_path = outdir / "validate_airr.tsv"
        validation = scrfu.tl.validate_airr(adata, airr_key=args.airr_key)
        validation.to_csv(validate_path, sep="\t", index=False)
        outputs.append(validate_path)

        manifest_path = outdir / "run_manifest.json"

        if args.dry_run:
            notes.append("Dry run: RFU calling, summaries, plots, and annotated h5ad were skipped.")
            outputs.append(manifest_path)
            _write_manifest(
                path=manifest_path,
                scrfu=scrfu,
                input_path=input_path,
                outdir=outdir,
                rfu_dir=args.rfu_dir,
                rfu_run=False,
                groupby=args.groupby,
                cell_type_key=args.cell_type_key,
                n_cells=adata.n_obs,
                output_files=outputs,
                notes=notes,
            )
            return 0

        _require_obs_columns(adata, [args.groupby])

        rfu_run = False
        if args.skip_rfu:
            _require_obs_columns(adata, ["rfu_label", "rfu_score"])
            notes.append("RFU execution skipped; existing adata.obs RFU columns were used.")
        else:
            scrfu.tl.call_rfu(
                adata,
                backend="rfu_repo",
                rfu_dir=args.rfu_dir,
                airr_key=args.airr_key,
                workdir=outdir / "rfu_work",
            )
            rfu_run = True

        global_summary_path = outdir / "rfu_summary_global.tsv"
        group_bar_path = outdir / "rfu_bar_by_group.png"
        score_hist_path = outdir / "rfu_score_hist.png"

        scrfu.tl.rfu_summary(adata).to_csv(global_summary_path, sep="\t", index=False)
        outputs.append(global_summary_path)

        outputs.extend(
            _write_group_outputs(
                adata,
                groupby=args.groupby,
                label="group",
                outdir=outdir,
                scrfu=scrfu,
                plt=plt,
            )
        )
        _save_plot(group_bar_path, lambda: scrfu.pl.rfu_bar(adata, groupby=args.groupby), plt)
        _save_plot(score_hist_path, lambda: scrfu.pl.rfu_score_hist(adata), plt)
        outputs.extend([group_bar_path, score_hist_path])

        if args.cell_type_key:
            if args.cell_type_key in adata.obs.columns:
                outputs.extend(
                    _write_group_outputs(
                        adata,
                        groupby=args.cell_type_key,
                        label="cell_type",
                        outdir=outdir,
                        scrfu=scrfu,
                        plt=plt,
                    )
                )
            else:
                notes.append(
                    f"cell_type_key '{args.cell_type_key}' not found; cell-type outputs skipped."
                )

        if args.write_annotated:
            annotated_path = outdir / "annotated.h5ad"
            adata.write_h5ad(annotated_path)
            outputs.append(annotated_path)

        outputs.append(manifest_path)
        _write_manifest(
            path=manifest_path,
            scrfu=scrfu,
            input_path=input_path,
            outdir=outdir,
            rfu_dir=args.rfu_dir,
            rfu_run=rfu_run,
            groupby=args.groupby,
            cell_type_key=args.cell_type_key,
            n_cells=adata.n_obs,
            output_files=outputs,
            notes=notes,
        )
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        return _fail(str(exc))

    print(f"Done. Outputs written to: {outdir}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run_benchmark(args)


if __name__ == "__main__":
    raise SystemExit(main())
