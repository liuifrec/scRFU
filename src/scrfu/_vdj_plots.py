"""Matplotlib-only drawing for the public plotting namespace."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from .visualization import (
    EmbeddingPlotData,
    LongitudinalPresenceResult,
    MembershipResult,
    SpectratypeResult,
    VJUsageResult,
    _bool,
    _ids,
    prepare_rfu_cdr3_spectratype,
    prepare_rfu_embedding,
    prepare_rfu_longitudinal_presence,
    prepare_rfu_membership,
    prepare_rfu_vj_usage,
)

if TYPE_CHECKING:
    from matplotlib.axes import Axes
else:
    Axes = Any

_COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#444444")


def _axes(ax: Axes | None) -> Axes:
    from .pl import _get_ax

    return _get_ax(ax)


def _palette(keys: Sequence[Any], palette: Mapping[Any, Any] | None) -> dict[Any, Any]:
    if palette is not None:
        missing = [k for k in keys if k not in palette]
        if missing:
            raise ValueError(f"Palette is missing categories: {missing}")
        return dict(palette)
    if len(keys) <= len(_COLORS):
        return dict(zip(keys, _COLORS, strict=False))
    from matplotlib import colormaps

    colors = colormaps["cividis"](np.linspace(0.04, 0.96, len(keys)))
    return dict(zip(keys, colors, strict=True))


def _matrix(
    ax: Axes,
    values: pd.DataFrame,
    *,
    cmap: str,
    colorbar: bool,
    value_label: str,
    vmin: float | None = None,
    vmax: float | None = None,
) -> Axes:
    if values.empty:
        raise ValueError("No observations are available for this plot.")
    image = ax.pcolormesh(
        np.arange(values.shape[1] + 1) - 0.5,
        np.arange(values.shape[0] + 1) - 0.5,
        values.to_numpy(dtype=float),
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        shading="flat",
        rasterized=False,
    )
    ax.set_xlim(-0.5, values.shape[1] - 0.5)
    ax.set_ylim(values.shape[0] - 0.5, -0.5)
    ax.set_xticks(range(values.shape[1]), [str(x) for x in values.columns], rotation=90)
    ax.set_yticks(range(values.shape[0]), [str(x) for x in values.index])
    if colorbar:
        bar = ax.figure.colorbar(image, ax=ax, label=value_label)
        bar.solids.set_rasterized(False)
    return ax


def _ribbon(
    ax: Axes,
    x0: float,
    x1: float,
    low0: float,
    high0: float,
    low1: float,
    high1: float,
    color: Any,
    alpha: float = 0.5,
) -> None:
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path

    dx = (x1 - x0) * 0.45
    vertices = [
        (x0, low0),
        (x0 + dx, low0),
        (x1 - dx, low1),
        (x1, low1),
        (x1, high1),
        (x1 - dx, high1),
        (x0 + dx, high0),
        (x0, high0),
        (x0, low0),
    ]
    codes = [
        Path.MOVETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
        Path.LINETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
        Path.CLOSEPOLY,
    ]
    ax.add_patch(PathPatch(Path(vertices, codes), facecolor=color, edgecolor="none", alpha=alpha))


def rfu_umap(
    data: Any = None,
    *,
    coordinates: Any = None,
    metadata: pd.DataFrame | None = None,
    cell_ids: Sequence[Any] | None = None,
    cell_col: str | None = None,
    basis: str = "X_umap",
    modality: str | None = None,
    allow_missing_metadata: bool = False,
    color: str = "rfu_label",
    highlight: Sequence[Any] | None = None,
    qualification_col: str | None = None,
    eligible_col: str | None = None,
    display_cells: Sequence[Any] | None = None,
    palette: Mapping[Any, Any] | None = None,
    labels: Mapping[Any, str] | None = None,
    background_color: Any = "#D8DBDF",
    size: float = 3,
    highlight_size: float | None = None,
    alpha: float = 1,
    legend: bool = True,
    title: str | None = None,
    ax: Axes | None = None,
) -> Axes:
    """Plot saved RNA/cellular coordinates, never an implicitly fitted embedding.

    Accepts AnnData, explicit MuData modality, indexed coordinates plus metadata,
    or ``tl.prepare_rfu_embedding`` output. RFUs or any categorical metadata column
    can provide color. An explicit qualification column instead displays saved
    boolean status; no score threshold is inferred. Missing values stay explicit.

    ``display_cells`` selects identifiers for display only, preserving coordinate
    order. It never determines a biological denominator. Selected highlights are
    overlaid on the same full display in original order. No file is saved or shown.
    """
    prepared = (
        data
        if isinstance(data, EmbeddingPlotData)
        else prepare_rfu_embedding(
            data,
            coordinates=coordinates,
            metadata=metadata,
            cell_ids=cell_ids,
            cell_col=cell_col,
            basis=basis,
            modality=modality,
            allow_missing_metadata=allow_missing_metadata,
        )
    )
    meta = prepared.metadata
    if qualification_col is not None:
        if qualification_col not in meta:
            raise ValueError(f"Missing qualification column {qualification_col!r}.")
        status = _bool(meta[qualification_col], qualification_col)
        values = status.map({True: "Qualified", False: "Not qualified"}).fillna(
            "Missing qualification"
        )
        if eligible_col is not None:
            if eligible_col not in meta:
                raise ValueError(f"Missing eligibility column {eligible_col!r}.")
            eligible = _bool(meta[eligible_col], eligible_col)
            if (status.fillna(False) & eligible.eq(False).fillna(False)).any():
                raise ValueError("A qualified assignment is marked ineligible.")
            values = values.mask(eligible.eq(False).fillna(False), "No eligible receptor")
    else:
        if color not in meta:
            raise ValueError(f"Embedding metadata is missing color column {color!r}.")
        if meta[color].dropna().eq("Missing").any() and meta[color].isna().any():
            raise ValueError("Literal color category 'Missing' collides with missing values.")
        values = meta[color].astype(object).mask(meta[color].isna(), "Missing")
    mask = np.ones(len(meta), dtype=bool)
    if display_cells is not None:
        selected = _ids(display_cells, context="Display cells")
        if not selected.isin(meta.index).all():
            raise ValueError("Display cell identifiers are absent from the original embedding.")
        mask = meta.index.isin(selected)
        if not mask.any():
            raise ValueError("The display selection is empty.")
    xy = prepared.coordinates[mask]
    values = values.loc[mask]
    keys = sorted(values.unique(), key=lambda x: (type(x).__name__, str(x)))
    if highlight is not None:
        keys = _ids(highlight, context="Highlight categories").tolist()
        if not set(keys).issubset(set(values)):
            raise ValueError("A highlight category has no cells in the display.")
    colors = _palette(keys, palette)
    ax = _axes(ax)
    if highlight is None:
        ax.scatter(
            xy[:, 0],
            xy[:, 1],
            c=[colors[v] for v in values],
            s=size,
            alpha=alpha,
            edgecolors="none",
        )
    else:
        ax.scatter(xy[:, 0], xy[:, 1], c=background_color, s=size, alpha=alpha, edgecolors="none")
        for key in keys:
            selected = values.eq(key).to_numpy()
            ax.scatter(
                xy[selected, 0],
                xy[selected, 1],
                c=[colors[key]],
                s=highlight_size or size,
                edgecolors="none",
            )
    if legend:
        from matplotlib.lines import Line2D

        handles = [
            Line2D(
                [],
                [],
                marker="o",
                linestyle="",
                color=colors[k],
                label=(labels or {}).get(k, str(k)),
            )
            for k in keys
        ]
        ax.legend(
            handles=handles, frameon=False, title="Qualification" if qualification_col else color
        )
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(title or "Existing embedding")
    return ax


def rfu_vj_usage(
    data: Any,
    *,
    groupby: str | None = None,
    unit: str | None = None,
    identity_cols: Sequence[str] | None = None,
    v_col: str = "v_call",
    j_col: str = "j_call",
    top_n: int | None = None,
    group_order: Sequence[Any] | None = None,
    segment: str = "v",
    style: str = "heatmap",
    group: Any = None,
    normalize: bool = True,
    cmap: str = "cividis",
    palette: Mapping[Any, Any] | None = None,
    labels: Mapping[Any, str] | None = None,
    colorbar: bool = True,
    vmin: float | None = None,
    vmax: float | None = None,
    title: str | None = None,
    ax: Axes | None = None,
) -> Axes:
    """Draw complete V/J marginals or a single group's joint V-to-J ribbons.

    ``segment='v'`` or ``'j'`` selects the heatmap marginal. A ribbon view uses
    both genes and requires an explicit ``group`` when multiple groups exist.
    Raw data follow ``tl.prepare_rfu_vj_usage``; prepared VJUsageResult is also
    accepted. Normalization divides by complete group weight, including missing
    and Other categories. No top-N limit is applied by default. Return the Axes.
    """
    if style not in {"heatmap", "ribbon"} or segment not in {"v", "j"}:
        raise ValueError("style must be 'heatmap'/'ribbon' and segment must be 'v'/'j'.")
    result = (
        data
        if isinstance(data, VJUsageResult)
        else prepare_rfu_vj_usage(
            data,
            groupby=groupby,
            unit=unit,
            identity_cols=identity_cols,
            v_col=v_col,
            j_col=j_col,
            top_n=top_n,
            group_order=group_order,
        )
    )
    if unit is not None and unit != result.unit:
        raise ValueError("unit disagrees with prepared data.")
    ax = _axes(ax)
    value_label = (
        f"Fraction of {result.unit} weight (complete group)"
        if normalize
        else f"{result.unit.title()} count"
    )
    if style == "heatmap":
        counts = result.v_counts if segment == "v" else result.j_counts
        values = (
            counts.div(result.denominators.replace(0, np.nan), axis=0)
            if normalize
            else counts.copy()
        )
        values.index = [
            (labels or {}).get(g, str(g)) + f" (n={result.denominators[g]:g})" for g in counts.index
        ]
        _matrix(
            ax, values, cmap=cmap, colorbar=colorbar, value_label=value_label, vmin=vmin, vmax=vmax
        )
        ax.set_xlabel(f"Literal {segment.upper()} call (missing and Other retained)")
        ax.set_ylabel("Group")
    else:
        from matplotlib.patches import Rectangle

        if group is None:
            if len(result.denominators) != 1:
                raise ValueError("Choose group explicitly for a multi-group ribbon plot.")
            group = result.denominators.index[0]
        if group not in result.denominators.index:
            raise ValueError("Ribbon group is not present in prepared data.")
        den = float(result.denominators[group])
        if den <= 0:
            raise ValueError("Ribbon group has zero observation weight.")
        scale = den if normalize else 1
        v, j = result.v_counts.loc[group] / scale, result.j_counts.loc[group] / scale
        colors = _palette(list(v.index), palette)
        v_offset = (v.cumsum() - v).to_dict()
        j_offset = (j.cumsum() - j).to_dict()
        cursors_v, cursors_j = v_offset.copy(), j_offset.copy()
        pairs = result.pairs.loc[result.pairs.group.eq(group)].set_index(["v_call", "j_call"]).value
        for vc in v.index:
            for jc in j.index:
                amount = float(pairs.get((vc, jc), 0)) / scale
                if amount > 0:
                    _ribbon(
                        ax,
                        0,
                        1,
                        cursors_v[vc],
                        cursors_v[vc] + amount,
                        cursors_j[jc],
                        cursors_j[jc] + amount,
                        colors[vc],
                        0.45,
                    )
                    cursors_v[vc] += amount
                    cursors_j[jc] += amount
        for side, counts, offsets in [(0, v, v_offset), (1, j, j_offset)]:
            for key, amount in counts.items():
                if amount > 0:
                    ax.add_patch(
                        Rectangle(
                            (side - 0.025, offsets[key]),
                            0.05,
                            amount,
                            facecolor=colors[key] if side == 0 else "#666666",
                            edgecolor="white",
                            linewidth=0.5,
                        )
                    )
                    ax.text(
                        side + (-0.05 if side == 0 else 0.05),
                        offsets[key] + amount / 2,
                        str(key),
                        ha="right" if side == 0 else "left",
                        va="center",
                    )
        ax.set_xlim(-0.35, 1.35)
        ax.set_ylim(0, den / scale)
        ax.set_xticks([0, 1], ["V call", "J call"])
        ax.set_yticks([])
        ax.set_ylabel(value_label)
        ax.set_title(title or f"{group}: V–J usage (n={den:g} {result.unit}s)")
        return ax
    ax.set_title(title or f"{segment.upper()} usage")
    return ax


def rfu_cdr3_spectratype(
    data: Any,
    *,
    groupby: str | None = None,
    unit: str | None = None,
    identity_cols: Sequence[str] | None = None,
    sequence_col: str = "cdr3aa",
    bins: int | Sequence[float] | None = None,
    invalid: str = "error",
    group_order: Sequence[Any] | None = None,
    normalize: bool = True,
    palette: Mapping[Any, Any] | None = None,
    labels: Mapping[Any, str] | None = None,
    legend: bool = True,
    linestyle: str = "-",
    linewidth: float = 1.5,
    title: str | None = None,
    ax: Axes | None = None,
) -> Axes:
    """Draw actual amino-acid CDR3 length distributions, with explicit valid denominators.

    Accepts raw/default-prepared data or SpectratypeResult. Missing, invalid and
    duplicate handling is defined by ``tl.prepare_rfu_cdr3_spectratype``. Labels
    report valid-sequence weight and excluded weight. Bins cannot hide tails.
    """
    result = (
        data
        if isinstance(data, SpectratypeResult)
        else prepare_rfu_cdr3_spectratype(
            data,
            groupby=groupby,
            unit=unit,
            identity_cols=identity_cols,
            sequence_col=sequence_col,
            bins=bins,
            invalid=invalid,
            group_order=group_order,
        )
    )
    if unit is not None and unit != result.unit:
        raise ValueError("unit disagrees with prepared data.")
    colors = _palette(list(result.counts.index), palette)
    ax = _axes(ax)
    for key, row in result.counts.iterrows():
        den = result.denominators[key]
        values = (
            row.to_numpy(dtype=float) / (den if den else np.nan)
            if normalize
            else row.to_numpy(dtype=float)
        )
        label = (labels or {}).get(
            key, str(key)
        ) + f" (n={den:g}; excluded={result.excluded.loc[key].sum():g})"
        ax.stairs(
            values,
            result.bin_edges,
            color=colors[key],
            label=label,
            linestyle=linestyle,
            linewidth=linewidth,
        )
    ax.set_xlabel("CDR3 length (amino acids)")
    ax.set_ylabel(
        f"Fraction of valid {result.unit} weight" if normalize else f"{result.unit.title()} count"
    )
    ax.set_title(title or "Observed CDR3 spectratype")
    if legend:
        ax.legend(frameon=False)
    return ax


def _membership_selection(
    result: MembershipResult, max_per_rfu: int | None, receptor_nodes: Sequence[str] | None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if max_per_rfu is not None and receptor_nodes is not None:
        raise ValueError("Choose max_per_rfu or receptor_nodes, not both.")
    nodes = result.receptors
    if receptor_nodes is not None:
        selected = _ids(receptor_nodes, context="Displayed receptor nodes")
        if not selected.isin(nodes.node_id).all():
            raise ValueError("Unknown receptor node in display selection.")
        nodes = nodes.set_index("node_id").loc[selected.tolist()].reset_index()
    elif max_per_rfu is not None:
        if isinstance(max_per_rfu, bool) or not isinstance(max_per_rfu, int) or max_per_rfu < 1:
            raise ValueError("max_per_rfu must be a positive integer or None.")
        chosen = set()
        for key in result.rfus.node_id:
            ids = result.edges.loc[result.edges.rfu_node.eq(key), "receptor_node"]
            z = nodes.loc[nodes.node_id.isin(ids)].sort_values(
                ["abundance", "node_id"], ascending=[False, True], kind="stable"
            )
            chosen.update(z.head(max_per_rfu).node_id)
        nodes = nodes.loc[nodes.node_id.isin(chosen)]
    if nodes.empty:
        raise ValueError("No receptor nodes selected for display.")
    return nodes, result.edges.loc[result.edges.receptor_node.isin(nodes.node_id)].copy()


def rfu_membership_network(
    data: Any,
    *,
    identity_cols: Sequence[str] = ("cdr3aa", "v_call"),
    rfu_col: str = "rfu_label",
    unit: str | None = None,
    donor_col: str | None = None,
    allow_multiple_rfus: bool = False,
    layout: str = "bipartite",
    max_per_rfu: int | None = None,
    receptor_nodes: Sequence[str] | None = None,
    positions: Mapping[str, Sequence[float]] | None = None,
    palette: Mapping[Any, Any] | None = None,
    node_color: Any = "#4A6580",
    node_colors: Mapping[str, Any] | None = None,
    open_singletons: bool = False,
    size_scale: float | None = None,
    rfu_size: float = 100,
    labels: Mapping[Any, str] | None = None,
    receptor_labels: bool = False,
    legend: bool = True,
    title: str | None = None,
    ax: Axes | None = None,
) -> Axes:
    """Draw confirmed observed-identity/RFU membership, not sequence-distance edges.

    Circle area is ``size_scale * abundance`` in the prepared count unit. When
    omitted, the scale caps the largest circle at 300 pt² (at most 12 pt²/unit).
    RFU squares have fixed area. Display limits never alter complete result
    counts. Radial layout requires one RFU per identity; bipartite layout permits
    explicitly allowed multi-membership. Optional positions use prepared node IDs.
    Neither geometry nor edges imply lineage, shared function or antigen binding.
    """
    if layout not in {"bipartite", "radial"}:
        raise ValueError("layout must be 'bipartite' or 'radial'.")
    result = (
        data
        if isinstance(data, MembershipResult)
        else prepare_rfu_membership(
            data,
            identity_cols=identity_cols,
            rfu_col=rfu_col,
            unit=unit,
            donor_col=donor_col,
            allow_multiple_rfus=allow_multiple_rfus,
        )
    )
    if unit is not None and unit != result.unit:
        raise ValueError("unit disagrees with prepared data.")
    nodes, edges = _membership_selection(result, max_per_rfu, receptor_nodes)
    rfus = result.rfus.loc[result.rfus.node_id.isin(edges.rfu_node)]
    colors = _palette(list(rfus.rfu), palette)
    pos = {}
    if positions is not None:
        for key in [*nodes.node_id, *rfus.node_id]:
            try:
                coordinate = np.asarray(positions[key], dtype=float)
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(
                    "positions must provide finite 2-D coordinates for every displayed node."
                ) from exc
            if coordinate.shape != (2,) or not np.isfinite(coordinate).all():
                raise ValueError(
                    "positions must provide finite 2-D coordinates for every displayed node."
                )
            pos[key] = tuple(coordinate)
    elif layout == "bipartite":
        pos.update(
            {
                key: (0, y)
                for key, y in zip(nodes.node_id, np.linspace(1, 0, len(nodes)), strict=True)
            }
        )
        pos.update(
            {key: (1, y) for key, y in zip(rfus.node_id, np.linspace(1, 0, len(rfus)), strict=True)}
        )
    else:
        if edges.groupby("receptor_node").rfu_node.nunique().gt(1).any():
            raise ValueError("Radial layout requires unique membership; use bipartite layout.")
        for i, key in enumerate(rfus.node_id):
            pos[key] = (i * 2.5, 0)
            ids = nodes.loc[
                nodes.node_id.isin(edges.loc[edges.rfu_node.eq(key), "receptor_node"]), "node_id"
            ].tolist()
            for radius, indices in zip(
                [0.5, 0.76, 1], np.array_split(np.arange(len(ids)), 3), strict=True
            ):
                for j, k in enumerate(indices):
                    angle = 2 * np.pi * j / len(indices) + np.pi / 2
                    pos[ids[k]] = (i * 2.5 + radius * np.cos(angle), radius * np.sin(angle))
    scale = (
        min(12, 300 / max(float(nodes.abundance.max()), 1)) if size_scale is None else size_scale
    )
    if not np.isfinite(scale) or scale <= 0 or not np.isfinite(rfu_size) or rfu_size <= 0:
        raise ValueError("Node area scales must be finite and positive.")
    ax = _axes(ax)
    for edge in edges.itertuples():
        a, b = pos[edge.receptor_node], pos[edge.rfu_node]
        (line,) = ax.plot([a[0], b[0]], [a[1], b[1]], color="#AEB4BD", linewidth=0.65, zorder=1)
        line.set_gid(f"{edge.receptor_node}->{edge.rfu_node}")
    for node in nodes.itertuples():
        col = (node_colors or {}).get(node.node_id, node_color)
        x, y = pos[node.node_id]
        artist = ax.scatter(
            [x],
            [y],
            s=scale * node.abundance,
            facecolors="white" if open_singletons and node.abundance == 1 else col,
            edgecolors=col,
            linewidths=0.5,
            zorder=3,
        )
        artist.set_gid(node.node_id)
        if receptor_labels:
            ax.annotate(node.node_id, (x, y), xytext=(-3, 3), textcoords="offset points")
    for node in rfus.itertuples():
        x, y = pos[node.node_id]
        artist = ax.scatter([x], [y], marker="s", s=rfu_size, color=colors[node.rfu], zorder=4)
        artist.set_gid(node.node_id)
        ax.annotate(
            (labels or {}).get(node.rfu, str(node.rfu)),
            (x, y),
            xytext=(8, 0),
            textcoords="offset points",
            va="center",
            color=colors[node.rfu],
        )
    ax.margins(0.15)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_axis_off()
    ax.set_title(title or "Observed identity → RFU membership")
    if legend:
        ax.text(
            0,
            -0.04,
            f"Circles: {scale:g} pt²/{result.unit}; squares: RFUs ({rfu_size:g} pt²).\nDisplayed identities: {len(nodes)}/{len(result.receptors)}; edges: saved membership only.",
            transform=ax.transAxes,
            va="top",
            fontsize="small",
        )
    return ax


def rfu_longitudinal_presence(
    data: Any,
    *,
    sample_col: str | None = None,
    sample_order: Sequence[Any] | None = None,
    identity_cols: Sequence[str] = ("rfu_label",),
    donor_col: str | None = None,
    unit: str | None = None,
    min_count: float = 1,
    sample_metadata: pd.DataFrame | None = None,
    missing_identity: str = "exclude",
    donor: Any = None,
    style: str = "heatmap",
    measure: str = "presence",
    normalize: bool = False,
    max_identities: int | None = None,
    cmap: str = "cividis",
    palette: Mapping[Any, Any] | None = None,
    colorbar: bool = True,
    labels: Mapping[Any, str] | None = None,
    title: str | None = None,
    ax: Axes | None = None,
) -> Axes:
    """Plot donor-stratified identity presence/abundance without implying clonal persistence.

    Heatmaps support any explicit sample order. Presence is detection at min_count;
    abundance optionally divides by included observation weight per sample. Missing
    visits stay blank. A two-sample ``style='alluvial'`` view counts distinct
    detected identities in the union and is available only for presence. It never
    tracks cells or weights ribbons by reads. Multi-donor inputs require an explicit
    donor selection; call repeatedly on supplied axes to compose donor panels.
    """
    if style not in {"heatmap", "alluvial"} or measure not in {"presence", "abundance"}:
        raise ValueError("Choose heatmap/alluvial style and presence/abundance measure.")
    if measure == "presence" and normalize:
        raise ValueError("Presence is binary; normalization applies only to abundance heatmaps.")
    if style == "alluvial" and (measure != "presence" or max_identities is not None):
        raise ValueError(
            "Alluvials require complete presence unions without identity downselection."
        )
    if isinstance(data, LongitudinalPresenceResult):
        result = data
    else:
        if sample_col is None or sample_order is None:
            raise ValueError("Raw longitudinal data require sample_col and sample_order.")
        result = prepare_rfu_longitudinal_presence(
            data,
            sample_col=sample_col,
            sample_order=sample_order,
            identity_cols=identity_cols,
            donor_col=donor_col,
            unit=unit,
            min_count=min_count,
            sample_metadata=sample_metadata,
            missing_identity=missing_identity,
        )
    if unit is not None and unit != result.unit:
        raise ValueError("unit disagrees with prepared data.")
    if donor is None:
        if len(result.matrices) != 1:
            raise ValueError("Select donor explicitly; donors are not silently pooled.")
        donor = next(iter(result.matrices))
    if donor not in result.matrices:
        raise ValueError("Unknown donor selection.")
    matrix = result.matrices[donor].copy()
    ax = _axes(ax)
    if style == "alluvial":
        from matplotlib.patches import Rectangle

        if len(result.sample_order) != 2 or not result.observed.loc[donor].all():
            raise ValueError(
                "Alluvial presence requires exactly two observed samples; missing visits are not zeros."
            )
        flows = result.presence_flows.loc[result.presence_flows.donor.eq(donor)].set_index("status")
        den = float(flows.denominator.iloc[0])
        if den <= 0:
            raise ValueError("No detected identities are available in this paired union.")
        statuses = ["First only", "Both", "Last only"]
        colors = _palette(
            statuses,
            palette or dict(zip(statuses, [_COLORS[1], _COLORS[2], _COLORS[0]], strict=True)),
        )
        a, b, c = [float(flows.loc[s, "identities"]) / den for s in statuses]
        for status, l0, h0, l1, h1 in [
            (statuses[0], 0, a, b + c, 1),
            (statuses[1], a, a + b, 0, b),
            (statuses[2], a + b, 1, b, b + c),
        ]:
            _ribbon(ax, 0, 1, l0, h0, l1, h1, colors[status], 0.6)
        for side, top in [(0, a + b), (1, b + c)]:
            ax.add_patch(Rectangle((side - 0.025, 0), 0.05, top, facecolor="#4F5967"))
            ax.add_patch(
                Rectangle(
                    (side - 0.025, top),
                    0.05,
                    1 - top,
                    facecolor="#DDDFE2",
                    hatch="///",
                    edgecolor="white",
                )
            )
        from matplotlib.lines import Line2D

        ax.legend(
            handles=[
                Line2D(
                    [], [], lw=5, color=colors[s], label=f"{s}: {int(flows.loc[s, 'identities'])}"
                )
                for s in statuses
            ],
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.12),
            ncol=3,
        )
        ax.set_xlim(-0.06, 1.06)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xticks([0, 1], [(labels or {}).get(s, str(s)) for s in result.sample_order])
        ax.set_ylabel(f"Detected identity union (n={den:g}); hatched = not observed")
    else:
        if max_identities is not None:
            if (
                isinstance(max_identities, bool)
                or not isinstance(max_identities, int)
                or max_identities < 1
            ):
                raise ValueError("max_identities must be a positive integer or None.")
            totals = matrix.sum(axis=1)
            chosen = sorted(
                range(len(matrix)), key=lambda i: (-totals.iloc[i], str(matrix.index[i]))
            )[:max_identities]
            matrix = matrix.iloc[chosen]
        if measure == "presence":
            values = matrix.ge(result.min_count).astype(float).mask(matrix.isna())
            value_label = f"Observed presence (≥{result.min_count:g} {result.unit})"
        else:
            values = (
                matrix.div(result.denominators.loc[donor].replace(0, np.nan), axis=1)
                if normalize
                else matrix
            )
            value_label = (
                f"Fraction of included {result.unit} weight"
                if normalize
                else f"{result.unit.title()} count"
            )
        values = values.copy()
        values.index = [f"Identity {i + 1}" for i in range(len(values))]
        values.columns = [
            (labels or {}).get(s, str(s))
            + (
                f" (n={result.denominators.loc[donor, s]:g})"
                if result.observed.loc[donor, s]
                else " (unavailable)"
            )
            for s in result.sample_order
        ]
        _matrix(
            ax,
            values,
            cmap=cmap,
            colorbar=colorbar,
            value_label=value_label,
            vmin=0,
            vmax=1 if measure == "presence" or normalize else None,
        )
        ax.set_xlabel("Sample/condition (explicit order)")
        ax.set_ylabel("Observed identity; display selection does not change denominators")
    ax.set_title(title or f"{donor}: observed identity {measure}")
    return ax
