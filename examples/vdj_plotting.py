"""Synthetic V(D)J/RFU gallery. No biological data or RFU assignment is used.

Run: python examples/vdj_plotting.py --outdir /tmp/scrfu-vdj-gallery
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import scrfu


def synthetic_receptors() -> pd.DataFrame:
    """Return twelve invented cells with known, deliberately mock memberships."""
    amino = [
        "CAAF",
        "CAAF",
        "CGGAF",
        "CAAF",
        "CTTQF",
        "CQQAAF",
        "CAAF",
        "CGGAF",
        "CQQAAF",
        "CAAF",
        "CTTQF",
        "CRRLGF",
    ]
    definitions = {
        "CAAF": ("TRBV1", "TRBJ1-1", "TGTGCTGCTTTT", "mock_rfu_A"),
        "CGGAF": ("TRBV2", "TRBJ1-1", "TGTGGTGGTGCTTTT", "mock_rfu_A"),
        "CTTQF": ("TRBV3", "TRBJ2-1", "TGTACTACTCAGTTT", "mock_rfu_B"),
        "CQQAAF": ("TRBV4", "TRBJ2-3", "TGTCAGCAGGCTGCTTTT", "mock_rfu_B"),
        "CRRLGF": ("TRBV2", "TRBJ2-1", "TGTCGTCGTCTGGGTTTT", "mock_rfu_B"),
    }
    frame = pd.DataFrame(
        {
            "input_row_id": [f"invented_record_{i:02d}" for i in range(12)],
            "cell_id": [f"invented_cell_{i:02d}" for i in range(12)],
            "chain": "TRB",
            "productive": True,
            "cdr3aa": amino,
            "v_call": [definitions[s][0] for s in amino],
            "j_call": [definitions[s][1] for s in amino],
            "junction": [definitions[s][2] for s in amino],
            "rfu_label": [definitions[s][3] for s in amino],
            "donor": ["invented_donor_a"] * 6 + ["invented_donor_b"] * 6,
            "visit": ["before"] * 3 + ["after"] * 3 + ["before"] * 3 + ["after"] * 3,
            "read_count": [4, 2, 5, 1, 3, 2, 8, 1, 1, 4, 2, 1],
            "pass_thr": True,
        }
    )
    frame.loc[1, "junction"] = "TGCGCCGCTTTC"
    return frame


def render_gallery(outdir: Path) -> None:
    """Compose supplied Axes and explicitly export the synthetic example."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = synthetic_receptors()
    # These are invented display coordinates, NOT a computed biological UMAP.
    coordinates = pd.DataFrame(
        [
            [-2, -1],
            [-1.7, -0.6],
            [-1.1, -0.9],
            [-1.4, 0.1],
            [0.7, 0.8],
            [1.3, 1.1],
            [-1.8, 0.5],
            [-0.9, 0.8],
            [1.7, 0.5],
            [-0.7, 0],
            [1.5, -0.3],
            [2, -0.6],
        ],
        index=rows.cell_id,
        columns=["x", "y"],
    )
    cell_data = scrfu.tl.prepare_receptor_plot_data(rows, unit="cell")
    vj = scrfu.tl.prepare_rfu_vj_usage(cell_data, groupby="rfu_label")
    network = scrfu.tl.prepare_rfu_membership(cell_data)
    presence = scrfu.tl.prepare_rfu_longitudinal_presence(
        cell_data,
        sample_col="visit",
        sample_order=["before", "after"],
        donor_col="donor",
        identity_cols=["cdr3aa", "v_call"],
    )
    with plt.rc_context(
        {
            "font.size": 8,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "svg.hashsalt": "scrfu-synthetic-vdj-gallery",
        }
    ):
        fig, axes = plt.subplots(2, 3, figsize=(14, 9), layout="constrained")
        scrfu.pl.rfu_umap(
            coordinates=coordinates,
            metadata=rows,
            cell_col="cell_id",
            size=35,
            ax=axes[0, 0],
            title="Invented coordinates; no biological embedding",
        )
        scrfu.pl.rfu_vj_usage(vj, segment="v", ax=axes[0, 1])
        scrfu.pl.rfu_vj_usage(
            cell_data, style="ribbon", ax=axes[0, 2], title="Synthetic joint V–J calls"
        )
        scrfu.pl.rfu_cdr3_spectratype(cell_data, groupby="rfu_label", ax=axes[1, 0])
        scrfu.pl.rfu_membership_network(network, open_singletons=True, ax=axes[1, 1])
        scrfu.pl.rfu_longitudinal_presence(
            presence, donor="invented_donor_a", style="alluvial", ax=axes[1, 2]
        )
        fig.suptitle(
            "scRFU visualization gallery — synthetic cells and mock RFU labels", fontsize=14
        )
        outdir.mkdir(parents=True, exist_ok=True)
        for suffix in ("svg", "pdf", "png"):
            metadata = (
                {"Date": None}
                if suffix == "svg"
                else {"CreationDate": None, "ModDate": None}
                if suffix == "pdf"
                else None
            )
            fig.savefig(outdir / f"synthetic_vdj_gallery.{suffix}", dpi=180, metadata=metadata)
        svg = outdir / "synthetic_vdj_gallery.svg"
        svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
        plt.close(fig)
    assert len(rows) == 12 and len(network.receptors) == 5
    assert np.array_equal(vj.denominators.to_numpy(), [7, 5])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    render_gallery(parser.parse_args().outdir)
