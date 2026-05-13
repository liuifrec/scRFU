from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes

from .summary import aggregate_rfu


def _get_ax(ax: Axes | None) -> Axes:
    if ax is not None:
        return ax
    _, ax = plt.subplots()
    return ax


def _select_top_columns(table: pd.DataFrame, top_n: int) -> pd.DataFrame:
    if table.empty or table.shape[1] == 0:
        return table

    top_n = max(int(top_n), 1)
    top = table.sum(axis=0).sort_values(ascending=False).head(top_n).index
    return table.loc[:, list(top)]


def _draw_empty(ax: Axes, title: str, ylabel: str) -> Axes:
    ax.set_title(title)
    ax.set_xlabel("RFU")
    ax.set_ylabel(ylabel)
    ax.text(0.5, 0.5, "No assigned RFU labels", transform=ax.transAxes, ha="center", va="center")
    return ax


def rfu_bar(
    adata: Any,
    groupby: str,
    top_n: int = 20,
    normalize: bool = True,
    ax: Axes | None = None,
) -> Axes:
    """
    Plot a grouped bar chart of the most abundant RFUs across groups.
    """
    ax = _get_ax(ax)
    table = aggregate_rfu(adata, groupby=groupby, normalize=normalize)
    table = _select_top_columns(table, top_n=top_n)

    ylabel = "Proportion" if normalize else "Count"
    if table.empty or table.shape[1] == 0:
        return _draw_empty(ax, "RFU abundance by group", ylabel)

    table.plot(kind="bar", ax=ax)
    ax.set_title("RFU abundance by group")
    ax.set_xlabel(groupby)
    ax.set_ylabel(ylabel)
    ax.legend(title="RFU", frameon=False)
    return ax


def rfu_heatmap(
    adata: Any,
    groupby: str,
    top_n: int = 50,
    normalize: bool = True,
    ax: Axes | None = None,
) -> Axes:
    """
    Plot a group-by-RFU heatmap for the most abundant RFUs.
    """
    ax = _get_ax(ax)
    table = aggregate_rfu(adata, groupby=groupby, normalize=normalize)
    table = _select_top_columns(table, top_n=top_n)

    ylabel = groupby
    if table.empty or table.shape[1] == 0:
        return _draw_empty(ax, "RFU heatmap", ylabel)

    im = ax.imshow(table.to_numpy(dtype=float), aspect="auto", interpolation="nearest")
    ax.set_title("RFU heatmap")
    ax.set_xlabel("RFU")
    ax.set_ylabel(groupby)
    ax.set_xticks(np.arange(table.shape[1]))
    ax.set_xticklabels(table.columns, rotation=90)
    ax.set_yticks(np.arange(table.shape[0]))
    ax.set_yticklabels(table.index.astype(str))
    ax.figure.colorbar(im, ax=ax, label="Proportion" if normalize else "Count")
    return ax


def rfu_score_hist(adata: Any, bins: int = 50, ax: Axes | None = None) -> Axes:
    """
    Plot a histogram of non-missing RFU scores.
    """
    if not hasattr(adata, "obs"):
        raise ValueError("AnnData-like object must have an .obs attribute.")
    if "rfu_score" not in adata.obs.columns:
        raise ValueError("adata.obs is missing required columns: ['rfu_score']")

    ax = _get_ax(ax)
    scores = pd.to_numeric(adata.obs["rfu_score"], errors="coerce").dropna()
    if scores.empty:
        ax.set_title("RFU score distribution")
        ax.set_xlabel("RFU score")
        ax.set_ylabel("Count")
        ax.text(
            0.5, 0.5, "No RFU scores available", transform=ax.transAxes, ha="center", va="center"
        )
        return ax

    ax.hist(scores.to_numpy(), bins=bins)
    ax.set_title("RFU score distribution")
    ax.set_xlabel("RFU score")
    ax.set_ylabel("Count")
    return ax


__all__ = ["rfu_bar", "rfu_heatmap", "rfu_score_hist"]
