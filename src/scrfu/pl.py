from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from .downstream import RFUOverlapResult, RFUPseudobulkResult
from .summary import aggregate_rfu
from .vdjdb import AntigenPermutationResult

if TYPE_CHECKING:
    from matplotlib.axes import Axes
else:
    Axes = Any


def _get_ax(ax: Axes | None) -> Axes:
    if ax is not None:
        return ax
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError(
            "Plotting requires matplotlib; install scRFU with `pip install scrfu[plotting]`."
        ) from exc
    _, ax = plt.subplots()
    return ax


def _heatmap(
    table: pd.DataFrame,
    *,
    ax: Axes,
    title: str,
    xlabel: str,
    ylabel: str,
    colorbar_label: str,
) -> Axes:
    image = ax.imshow(table.to_numpy(dtype=float), aspect="auto", interpolation="nearest")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xticks(np.arange(table.shape[1]))
    ax.set_xticklabels(table.columns.astype(str), rotation=90)
    ax.set_yticks(np.arange(table.shape[0]))
    ax.set_yticklabels(table.index.astype(str))
    ax.figure.colorbar(image, ax=ax, label=colorbar_label)
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


def rfu_metric_heatmap(
    data: RFUPseudobulkResult | pd.DataFrame,
    *,
    top_n: int | None = 50,
    row_standardize: bool = False,
    value_label: str | None = None,
    ax: Axes | None = None,
) -> Axes:
    """Plot RFUs by samples from a pseudobulk result or sample-by-RFU matrix."""
    matrix = data.matrix.copy() if isinstance(data, RFUPseudobulkResult) else data.copy()
    if matrix.empty or matrix.shape[1] == 0:
        raise ValueError("RFU metric heatmap input is empty.")
    matrix = matrix.apply(pd.to_numeric, errors="raise")
    if top_n is not None:
        if top_n < 1:
            raise ValueError("top_n must be positive or None.")
        totals = matrix.abs().sum(axis=0)
        selected = sorted(matrix.columns, key=lambda item: (-totals[item], str(item)))[:top_n]
        matrix = matrix.loc[:, selected]
    table = matrix.T
    if row_standardize:
        means = table.mean(axis=1)
        deviations = table.std(axis=1, ddof=0).replace(0, np.nan)
        table = table.sub(means, axis=0).div(deviations, axis=0).fillna(0.0)
    return _heatmap(
        table,
        ax=_get_ax(ax),
        title="RFU metric heatmap",
        xlabel=str(matrix.index.name or "sample"),
        ylabel=str(matrix.columns.name or "RFU"),
        colorbar_label=value_label or ("Row z-score" if row_standardize else "Value"),
    )


def rfu_overlap_heatmap(
    data: RFUOverlapResult | pd.DataFrame,
    *,
    metric: str | None = None,
    direction: str | None = None,
    ax: Axes | None = None,
) -> Axes:
    """Plot a square RFU overlap matrix with explicit similarity/distance labeling."""
    if isinstance(data, RFUOverlapResult):
        matrix = data.matrix.copy()
        metric = data.metric
        direction = data.direction
    else:
        matrix = data.copy()
    if matrix.empty:
        raise ValueError("RFU overlap heatmap input is empty.")
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("RFU overlap heatmap requires a square matrix.")
    if matrix.index.astype(str).tolist() != matrix.columns.astype(str).tolist():
        raise ValueError("RFU overlap heatmap row and column sample labels must match in order.")
    direction = direction or "similarity"
    return _heatmap(
        matrix,
        ax=_get_ax(ax),
        title=f"RFU overlap: {metric or 'metric'} ({direction})",
        xlabel="Sample",
        ylabel="Sample",
        colorbar_label=direction.title(),
    )


def rfu_convergence(
    metrics: pd.DataFrame,
    *,
    richness_col: str = "unique_cdr3_richness",
    abundance_col: str = "cell_abundance",
    rfu_col: str = "rfu_label",
    phenotype_col: str | None = None,
    annotate_top: int = 0,
    ax: Axes | None = None,
) -> Axes:
    """Plot RFU sequence richness against abundance or multiplicity."""
    required = [richness_col, abundance_col, rfu_col]
    if phenotype_col is not None:
        required.append(phenotype_col)
    missing = [column for column in required if column not in metrics]
    if missing:
        raise ValueError(f"RFU convergence input is missing columns: {missing}")
    if metrics.empty:
        raise ValueError("RFU convergence input is empty.")
    plot_ax = _get_ax(ax)
    if phenotype_col is None:
        plot_ax.scatter(metrics[richness_col], metrics[abundance_col])
    else:
        for phenotype in sorted(metrics[phenotype_col].dropna().unique(), key=str):
            subset = metrics.loc[metrics[phenotype_col].eq(phenotype)]
            plot_ax.scatter(subset[richness_col], subset[abundance_col], label=str(phenotype))
        plot_ax.legend(title=phenotype_col, frameon=False)
    if annotate_top:
        if annotate_top < 0:
            raise ValueError("annotate_top must be non-negative.")
        ranked = metrics.assign(_order=range(len(metrics))).sort_values(
            [abundance_col, richness_col, "_order"], ascending=[False, False, True], kind="stable"
        )
        for _, row in ranked.head(annotate_top).iterrows():
            plot_ax.annotate(str(row[rfu_col]), (row[richness_col], row[abundance_col]))
    plot_ax.set_xlabel(richness_col)
    plot_ax.set_ylabel(abundance_col)
    plot_ax.set_title("RFU convergence")
    return plot_ax


def rfu_phenotype_heatmap(
    coupling: pd.DataFrame,
    *,
    rfu_col: str = "rfu_label",
    phenotype_col: str = "phenotype",
    value_col: str = "phenotype_specific_proportion",
    top_n: int | None = 50,
    min_prevalence: float = 0.0,
    ax: Axes | None = None,
) -> Axes:
    """Plot RFU-by-phenotype coupling values from the long-form result table."""
    required = [rfu_col, phenotype_col, value_col]
    missing = [column for column in required if column not in coupling]
    if missing:
        raise ValueError(f"RFU phenotype heatmap input is missing columns: {missing}")
    if coupling.empty:
        raise ValueError("RFU phenotype heatmap input is empty.")
    table = coupling.pivot(index=rfu_col, columns=phenotype_col, values=value_col).fillna(0.0)
    prevalence = (table > 0).mean(axis=1)
    table = table.loc[prevalence >= min_prevalence]
    if top_n is not None:
        if top_n < 1:
            raise ValueError("top_n must be positive or None.")
        totals = table.abs().sum(axis=1)
        selected = sorted(table.index, key=lambda item: (-totals[item], str(item)))[:top_n]
        table = table.loc[selected]
    if table.empty:
        raise ValueError("No RFUs remain after phenotype heatmap filtering.")
    return _heatmap(
        table,
        ax=_get_ax(ax),
        title="RFU phenotype coupling",
        xlabel=phenotype_col,
        ylabel="RFU",
        colorbar_label=value_col,
    )


def repertoire_metric_comparison(
    repertoire: pd.DataFrame,
    rfu: pd.DataFrame,
    *,
    sample_key: str,
    repertoire_metric: str,
    rfu_metric: str,
    ax: Axes | None = None,
) -> Axes:
    """Scatter two descriptive sample-level metrics without inferential claims."""
    for name, table, metric_name in (
        ("repertoire", repertoire, repertoire_metric),
        ("RFU", rfu, rfu_metric),
    ):
        missing = [column for column in (sample_key, metric_name) if column not in table]
        if missing:
            raise ValueError(f"{name} metric input is missing columns: {missing}")
    merged = repertoire[[sample_key, repertoire_metric]].merge(
        rfu[[sample_key, rfu_metric]], on=sample_key, how="inner", validate="one_to_one"
    )
    if merged.empty:
        raise ValueError("Repertoire comparison has no matched samples.")
    plot_ax = _get_ax(ax)
    plot_ax.scatter(merged[repertoire_metric], merged[rfu_metric])
    for _, row in merged.iterrows():
        plot_ax.annotate(str(row[sample_key]), (row[repertoire_metric], row[rfu_metric]))
    plot_ax.set_xlabel(repertoire_metric)
    plot_ax.set_ylabel(rfu_metric)
    plot_ax.set_title("RFU and conventional repertoire metrics")
    return plot_ax


def _antigen_abundance_table(
    abundance: pd.DataFrame,
    *,
    rfu_col: str,
    antigen_col: str,
    value_col: str,
    top_n_rfus: int | None,
    top_n_antigens: int | None,
) -> pd.DataFrame:
    if abundance.empty:
        raise ValueError("RFU-antigen plot input is empty.")
    required = [rfu_col, antigen_col, value_col]
    missing = [column for column in required if column not in abundance]
    if missing:
        raise ValueError(f"RFU-antigen input is missing columns: {missing}")
    data = abundance.loc[
        abundance[rfu_col].notna() & abundance[antigen_col].notna(), required
    ].copy()
    if data.empty:
        raise ValueError("RFU-antigen plot input is empty.")
    data[value_col] = pd.to_numeric(data[value_col], errors="raise")
    totals = data.groupby(rfu_col, observed=True)[value_col].sum()
    antigens = data.groupby(antigen_col, observed=True)[value_col].sum()
    ordered_rfus = sorted(totals.index, key=lambda item: (-totals[item], str(item)))
    ordered_antigens = sorted(antigens.index, key=lambda item: (-antigens[item], str(item)))
    if top_n_rfus is not None:
        if top_n_rfus < 1:
            raise ValueError("top_n_rfus must be positive or None.")
        ordered_rfus = ordered_rfus[:top_n_rfus]
    if top_n_antigens is not None:
        if top_n_antigens < 1:
            raise ValueError("top_n_antigens must be positive or None.")
        ordered_antigens = ordered_antigens[:top_n_antigens]
    table = data.pivot_table(
        index=rfu_col,
        columns=antigen_col,
        values=value_col,
        aggfunc="sum",
        fill_value=0.0,
    )
    return table.reindex(index=ordered_rfus, columns=ordered_antigens, fill_value=0.0)


def rfu_antigen_heatmap(
    abundance: pd.DataFrame,
    *,
    rfu_col: str = "rfu_label",
    antigen_col: str = "epitope",
    value_col: str = "within_rfu_proportion",
    top_n_rfus: int | None = 30,
    top_n_antigens: int | None = 30,
    ax: Axes | None = None,
) -> Axes:
    """Plot deterministic RFU-by-antigen evidence abundance or proportion."""
    table = _antigen_abundance_table(
        abundance,
        rfu_col=rfu_col,
        antigen_col=antigen_col,
        value_col=value_col,
        top_n_rfus=top_n_rfus,
        top_n_antigens=top_n_antigens,
    )
    return _heatmap(
        table,
        ax=_get_ax(ax),
        title="RFU antigen evidence",
        xlabel=antigen_col,
        ylabel=rfu_col,
        colorbar_label=value_col,
    )


def rfu_antigen_coherence(
    metrics: pd.DataFrame,
    *,
    x: str = "vdjdb_matched_sequences",
    y: str = "antigen_purity",
    size: str = "total_rfu_sequences",
    color: str | None = None,
    ax: Axes | None = None,
) -> Axes:
    """Plot descriptive RFU antigen-evidence coherence without inferential claims."""
    required = [x, y, size, *([color] if color else [])]
    missing = [column for column in required if column not in metrics]
    if missing:
        raise ValueError(f"RFU antigen-coherence input is missing columns: {missing}")
    data = metrics.dropna(subset=[x, y, size]).copy()
    if data.empty:
        raise ValueError("RFU antigen-coherence plot has no finite rows.")
    x_values = pd.to_numeric(data[x], errors="raise")
    y_values = pd.to_numeric(data[y], errors="raise")
    sizes = pd.to_numeric(data[size], errors="raise").clip(lower=0)
    scaled_sizes = 25.0 + 125.0 * sizes / max(float(sizes.max()), 1.0)
    plot_ax = _get_ax(ax)
    if color is None:
        plot_ax.scatter(x_values, y_values, s=scaled_sizes, alpha=0.75)
    else:
        categorical = pd.Categorical(data[color].astype("string"))
        points = plot_ax.scatter(
            x_values,
            y_values,
            s=scaled_sizes,
            c=categorical.codes,
            alpha=0.75,
        )
        labels = [str(label) for label in categorical.categories]
        handles, _ = points.legend_elements()
        plot_ax.legend(handles, labels, title=color, frameon=False)
    plot_ax.set_xlabel(x)
    plot_ax.set_ylabel(y)
    plot_ax.set_title("RFU antigen-evidence coherence")
    return plot_ax


def antigen_permutation_distribution(
    result: AntigenPermutationResult,
    *,
    bins: int = 30,
    ax: Axes | None = None,
) -> Axes:
    """Plot an observed coherence statistic against its permutation null."""
    if not isinstance(result, AntigenPermutationResult):
        raise TypeError("result must be an AntigenPermutationResult.")
    values = np.asarray(result.permutation_values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        raise ValueError("Permutation result contains no finite null values.")
    if bins < 1:
        raise ValueError("bins must be positive.")
    plot_ax = _get_ax(ax)
    plot_ax.hist(values, bins=bins, alpha=0.7, label="Permutation null")
    plot_ax.axvline(result.observed, color="black", linestyle="--", label="Observed")
    plot_ax.set_xlabel(str(result.parameters.get("metric", "coherence statistic")))
    plot_ax.set_ylabel("Permutation count")
    plot_ax.set_title("RFU antigen-coherence null benchmark")
    plot_ax.legend(frameon=False)
    return plot_ax


def rfu_antigen_bubble(
    abundance: pd.DataFrame,
    *,
    rfu_col: str = "rfu_label",
    antigen_col: str = "epitope",
    size_col: str = "antigen_abundance",
    color_col: str = "within_rfu_proportion",
    top_n_rfus: int | None = 30,
    top_n_antigens: int | None = 30,
    ax: Axes | None = None,
) -> Axes:
    """Plot RFU-antigen evidence with abundance size and within-RFU proportion color."""
    size_table = _antigen_abundance_table(
        abundance,
        rfu_col=rfu_col,
        antigen_col=antigen_col,
        value_col=size_col,
        top_n_rfus=top_n_rfus,
        top_n_antigens=top_n_antigens,
    )
    color_table = _antigen_abundance_table(
        abundance,
        rfu_col=rfu_col,
        antigen_col=antigen_col,
        value_col=color_col,
        top_n_rfus=None,
        top_n_antigens=None,
    ).reindex(index=size_table.index, columns=size_table.columns, fill_value=0.0)
    sizes = size_table.to_numpy(dtype=float)
    scaled = 20.0 + 180.0 * sizes / max(float(sizes.max()), 1.0)
    y_positions, x_positions = np.indices(size_table.shape)
    plot_ax = _get_ax(ax)
    points = plot_ax.scatter(
        x_positions.ravel(),
        y_positions.ravel(),
        s=scaled.ravel(),
        c=color_table.to_numpy(dtype=float).ravel(),
    )
    plot_ax.set_xticks(np.arange(size_table.shape[1]))
    plot_ax.set_xticklabels(size_table.columns.astype(str), rotation=90)
    plot_ax.set_yticks(np.arange(size_table.shape[0]))
    plot_ax.set_yticklabels(size_table.index.astype(str))
    plot_ax.set_xlabel(antigen_col)
    plot_ax.set_ylabel(rfu_col)
    plot_ax.set_title("RFU antigen evidence")
    plot_ax.figure.colorbar(points, ax=plot_ax, label=color_col)
    return plot_ax


__all__ = [
    "antigen_permutation_distribution",
    "repertoire_metric_comparison",
    "rfu_antigen_bubble",
    "rfu_antigen_coherence",
    "rfu_antigen_heatmap",
    "rfu_bar",
    "rfu_convergence",
    "rfu_heatmap",
    "rfu_metric_heatmap",
    "rfu_overlap_heatmap",
    "rfu_phenotype_heatmap",
    "rfu_score_hist",
]
