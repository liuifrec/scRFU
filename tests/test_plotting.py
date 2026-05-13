from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import pytest
from matplotlib.axes import Axes

import scrfu


class DummyAdata:
    def __init__(self, obs: pd.DataFrame) -> None:
        self.obs = obs


def test_rfu_bar_returns_axes():
    adata = DummyAdata(
        pd.DataFrame(
            {
                "group": ["g1", "g1", "g2", "g2"],
                "rfu_label": ["RFU2", "RFU10", "RFU2", pd.NA],
                "rfu_score": [0.1, 0.2, 0.3, pd.NA],
            }
        )
    )

    ax = scrfu.pl.rfu_bar(adata, groupby="group", top_n=2, normalize=True)

    assert isinstance(ax, Axes)
    plt.close(ax.figure)


def test_rfu_heatmap_returns_axes():
    adata = DummyAdata(
        pd.DataFrame(
            {
                "group": ["g1", "g1", "g2"],
                "rfu_label": ["RFU1", "RFU2", "RFU10"],
                "rfu_score": [0.1, 0.2, 0.3],
            }
        )
    )

    ax = scrfu.pl.rfu_heatmap(adata, groupby="group", top_n=3, normalize=False)

    assert isinstance(ax, Axes)
    plt.close(ax.figure)


def test_rfu_score_hist_returns_axes():
    adata = DummyAdata(pd.DataFrame({"rfu_score": [0.1, pd.NA, 0.3, 0.5]}))

    ax = scrfu.pl.rfu_score_hist(adata, bins=5)

    assert isinstance(ax, Axes)
    plt.close(ax.figure)


def test_plotting_handles_all_unassigned_rfu_gracefully():
    adata = DummyAdata(
        pd.DataFrame(
            {
                "group": ["g1", "g2"],
                "rfu_label": [pd.NA, pd.NA],
                "rfu_score": [pd.NA, pd.NA],
            }
        )
    )

    bar_ax = scrfu.pl.rfu_bar(adata, groupby="group")
    heatmap_ax = scrfu.pl.rfu_heatmap(adata, groupby="group")
    hist_ax = scrfu.pl.rfu_score_hist(adata)

    assert isinstance(bar_ax, Axes)
    assert isinstance(heatmap_ax, Axes)
    assert isinstance(hist_ax, Axes)
    plt.close(bar_ax.figure)
    plt.close(heatmap_ax.figure)
    plt.close(hist_ax.figure)


def test_plotting_raises_for_missing_required_columns():
    adata = DummyAdata(pd.DataFrame({"group": ["g1"]}))

    with pytest.raises(ValueError, match="missing required columns"):
        scrfu.pl.rfu_bar(adata, groupby="group")

    with pytest.raises(ValueError, match="missing required columns"):
        scrfu.pl.rfu_heatmap(adata, groupby="group")

    with pytest.raises(ValueError, match="missing required columns"):
        scrfu.pl.rfu_score_hist(adata)
