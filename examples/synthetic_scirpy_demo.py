from __future__ import annotations

import os
import sys
from importlib import import_module
from pathlib import Path

import anndata as ad
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

DEFAULT_OUTDIR = REPO_ROOT / ".scrfu" / "demo"
DEMO_OUTDIR = Path(os.environ.get("SCRFU_DEMO_OUTDIR", DEFAULT_OUTDIR)).expanduser().resolve()
MPLCONFIGDIR = DEMO_OUTDIR / ".mplconfig"
MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))


def build_demo_adata() -> ad.AnnData:
    obs_index = [f"cell{i}" for i in range(1, 11)]
    obs = pd.DataFrame(
        {
            "sample": ["S1", "S1", "S1", "S1", "S1", "S2", "S2", "S2", "S2", "S2"],
            "cell_type": [
                "CD8_T",
                "CD4_T",
                "CD8_T",
                "Treg",
                "MAIT",
                "CD8_T",
                "CD4_T",
                "NK_like",
                "Treg",
                "GammaDelta_like",
            ],
        },
        index=pd.Index(obs_index),
    )
    adata = ad.AnnData(obs=obs)
    adata.obsm["airr"] = pd.DataFrame(
        {
            "cell_id": obs_index,
            "chain": ["TRB", "TRB", "TRB", "TRA", "TRB", "TRB", "TRB", "TRB", "TRB", "TRB"],
            "cdr3aa": [
                "CASSLGQETQYF",
                "CASSPPSGGNEQF",
                "CASSIRSSYEQYF",
                "CAVRDTNTNAGKSTF",
                "CASSQETQYF",
                "CASSLGQETQYF",
                "CASSDLAGNTIYF",
                "CASSPGLAGNTIYF",
                "CASSIRSSYEQYF",
                "CASSQGELFF",
            ],
            "v_call": [
                "TRBV7-9",
                "TRBV5-1",
                "TRBV6-5",
                "TRAV12-2",
                "TRBV20-1",
                "TRBV7-9",
                "TRBV4-1",
                "TRBV27",
                "TRBV6-5",
                "TRBV24-1",
            ],
            "productive": [True, True, True, True, False, True, True, True, True, False],
        },
        index=obs.index.copy(),
    )
    return adata


def build_synthetic_rfu(features: pd.DataFrame) -> pd.DataFrame:
    mapping = {
        "cell1": ("RFU10", 0.91),
        "cell2": ("RFU2", 0.73),
        "cell3": ("RFU10", 0.88),
        "cell6": ("RFU10", 0.95),
        "cell7": ("RFU3", 0.64),
        "cell8": ("RFU2", 0.69),
        "cell9": ("RFU10", 0.86),
    }
    rows = []
    for cell_id in features["cell_id"]:
        if cell_id in mapping:
            label, score = mapping[cell_id]
            rows.append({"cell_id": cell_id, "rfu_label": label, "rfu_score": score})
    return pd.DataFrame(rows)


def output_dir() -> Path:
    return DEMO_OUTDIR


def _import_runtime_modules() -> tuple[object, object, object, object, object, object]:
    matplotlib = import_module("matplotlib")
    matplotlib.use("Agg")
    plt = import_module("matplotlib.pyplot")
    attach = import_module("scrfu.attach")
    extract = import_module("scrfu.extract")
    io = import_module("scrfu.io")
    pl = import_module("scrfu.pl")
    tl = import_module("scrfu.tl")
    return plt, attach, extract, io, pl, tl


def save_plot(path: Path, plotter, plt: object) -> None:
    ax = plotter()
    ax.figure.tight_layout()
    ax.figure.savefig(path, dpi=150)
    plt.close(ax.figure)


def main() -> int:
    outdir = output_dir()
    outdir.mkdir(parents=True, exist_ok=True)
    plt, attach, extract, io, pl, tl = _import_runtime_modules()

    adata = build_demo_adata()
    features = extract.extract_trb_features(adata, airr_key="airr", chain="TRB")
    synthetic_rfu = build_synthetic_rfu(features)

    attach.attach_rfu_results(
        adata,
        features,
        synthetic_rfu,
        provenance={"backend": "synthetic_demo", "note": "No upstream RFU execution."},
    )

    global_summary = tl.rfu_summary(adata)
    grouped_summary = tl.rfu_summary(adata, groupby="sample")
    rfu_matrix = tl.aggregate_rfu(adata, groupby="sample")
    matrix_path = outdir / "rfu_matrix.tsv"
    io.export_rfu_matrix(adata, groupby="sample", path=matrix_path)

    bar_path = outdir / "rfu_bar.png"
    heatmap_path = outdir / "rfu_heatmap.png"
    hist_path = outdir / "rfu_score_hist.png"

    save_plot(bar_path, lambda: pl.rfu_bar(adata, groupby="sample"), plt)
    save_plot(heatmap_path, lambda: pl.rfu_heatmap(adata, groupby="sample"), plt)
    save_plot(hist_path, lambda: pl.rfu_score_hist(adata), plt)

    print("Global summary:")
    print(global_summary.to_string(index=False))
    print()
    print("Grouped summary:")
    print(grouped_summary.to_string(index=False))
    print()
    print("RFU matrix:")
    print(rfu_matrix.to_string())
    print()
    print("Wrote outputs:")
    for path in [matrix_path, bar_path, heatmap_path, hist_path]:
        print(path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
