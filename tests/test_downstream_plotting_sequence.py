from __future__ import annotations

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pytest
from matplotlib.axes import Axes

from scrfu.pl import (
    repertoire_metric_comparison,
    rfu_convergence,
    rfu_metric_heatmap,
    rfu_overlap_heatmap,
    rfu_phenotype_heatmap,
)
from scrfu.tl import rfu_overlap, rfu_sequence_matrix

matplotlib.use("Agg")


def test_downstream_plots_return_axes_and_labels(monkeypatch) -> None:
    monkeypatch.setattr(plt, "show", lambda: pytest.fail("plt.show() must not be called"))
    matrix = pd.DataFrame(
        [[2, 1], [0, 3]], index=pd.Index(["s1", "s2"], name="sample"), columns=["r1", "r2"]
    )
    ax1 = rfu_metric_heatmap(matrix, top_n=1)
    ax2 = rfu_overlap_heatmap(rfu_overlap(matrix, metric="jaccard"))
    metrics = pd.DataFrame(
        {"rfu_label": ["r1", "r2"], "unique_cdr3_richness": [2, 1], "cell_abundance": [4, 2]}
    )
    ax3 = rfu_convergence(metrics, annotate_top=1)
    coupling = pd.DataFrame(
        {
            "rfu_label": ["r1", "r1", "r2", "r2"],
            "phenotype": ["A", "B", "A", "B"],
            "phenotype_specific_proportion": [1.0, 0.0, 0.5, 0.5],
        }
    )
    ax4 = rfu_phenotype_heatmap(coupling)
    conventional = pd.DataFrame({"sample": ["s1", "s2"], "shannon": [0.2, 0.5]})
    rfu = pd.DataFrame({"sample": ["s1", "s2"], "richness": [2, 3]})
    ax5 = repertoire_metric_comparison(
        conventional, rfu, sample_key="sample", repertoire_metric="shannon", rfu_metric="richness"
    )
    assert all(isinstance(ax, Axes) for ax in (ax1, ax2, ax3, ax4, ax5))
    assert ax2.get_title().endswith("(similarity)")
    assert ax3.get_xlabel() == "unique_cdr3_richness"


@pytest.mark.parametrize(
    "function,args",
    [
        (rfu_metric_heatmap, (pd.DataFrame(),)),
        (rfu_overlap_heatmap, (pd.DataFrame(),)),
        (
            rfu_convergence,
            (pd.DataFrame(columns=["rfu_label", "unique_cdr3_richness", "cell_abundance"]),),
        ),
        (
            rfu_phenotype_heatmap,
            (pd.DataFrame(columns=["rfu_label", "phenotype", "phenotype_specific_proportion"]),),
        ),
    ],
)
def test_downstream_plots_reject_empty(function, args) -> None:
    with pytest.raises(ValueError, match="empty"):
        function(*args)


def test_sequence_matrix_alignments_and_weighting() -> None:
    frame = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3"],
            "cdr3aa": ["CAF", "CAF", "CASSW"],
            "rfu_label": ["r1", "r1", "r1"],
            "clonotype_id": ["a", "a", "b"],
        }
    )
    unique = rfu_sequence_matrix(
        frame, rfu="r1", weighting="unique_sequence", alignment="conserved_ends"
    )
    cell = rfu_sequence_matrix(frame, rfu="r1", weighting="cell", alignment="left")
    assert unique.shape == (5, 21)
    assert unique.loc[1, "C"] == 1
    assert unique.loc[5, "F"] == pytest.approx(0.5)
    assert unique.loc[5, "W"] == pytest.approx(0.5)
    assert cell.loc[2, "A"] == pytest.approx(1.0)
    with pytest.raises(ValueError, match="Invalid amino-acid"):
        rfu_sequence_matrix(pd.DataFrame({"cdr3aa": ["CA*F"]}))
