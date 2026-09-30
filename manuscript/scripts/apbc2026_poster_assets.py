"""Extract poster cutouts from immutable figure tables; no result calculations.

Run with --out /tmp/apbc-preview for review, or omit --out for the external pack.
Only this exporter and its asset manifest belong in Git. The optional completed
checkpoint in that manifest is excluded from the rendering identity to avoid a
self-referential hash. All source values are checked against rendered artists.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import shlex
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager, patches
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.ticker import PercentFormatter
from PIL import Image, ImageDraw, ImageFont

from ._rp1_14_figure_drawing import COLORS, GRAY, INK, layout_audit
from .rp1_14_reusable_figures import REPO, json_save, require, safe_child, sha256

MANIFEST = REPO / "manuscript/figures/apbc2026_poster_asset_manifest.json"
RFU = "#72549A"
STYLE = {
    "font.family": "DejaVu Sans",
    "font.size": 24,
    "text.color": INK,
    "axes.labelsize": 26,
    "axes.titlesize": 28,
    "axes.labelcolor": INK,
    "axes.edgecolor": GRAY,
    "axes.linewidth": 1.1,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.labelsize": 24,
    "ytick.labelsize": 24,
    "xtick.major.size": 5,
    "ytick.major.size": 0,
    "legend.fontsize": 24,
    "lines.markersize": 11,
    "lines.linewidth": 2,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "pdf.fonttype": 42,
    "svg.fonttype": "path",
    "svg.hashsalt": "apbc2026-poster-assets-v1",
    "figure.max_open_warning": 0,
}


def plain(value):
    """Normalize artist/source types for exact, NaN-safe comparison and JSON."""
    if isinstance(value, (list, tuple, np.ndarray)):
        return [plain(v) for v in value]
    if pd.isna(value):
        return None
    return value.item() if isinstance(value, np.generic) else value


class Trace:
    """Bind actual artist coordinates, colors and labels to frozen table cells."""

    def __init__(self, tables):
        self.tables = tables
        self.probes = []

    def bind(self, table, rows, columns, actual, formatter=None, role="coordinate"):
        self.probes.append((table, list(rows), list(columns), actual, formatter, role))

    def check(self):
        records = []
        for table, rows, columns, actual, formatter, role in self.probes:
            source = self.tables[table].loc[rows, columns]
            expected = formatter(source) if formatter else source.to_numpy().ravel()
            require(
                plain(actual()) == plain(expected),
                f"Artist/source mismatch: {table} rows={rows} columns={columns} role={role}",
            )
            records.append(
                {
                    "table": table,
                    "tsv_line_numbers": [int(i) + 2 for i in rows],
                    "columns": columns,
                    "frozen_values": plain(source.to_numpy()),
                    "rendered_values": plain(actual()),
                    "role": role,
                }
            )
        require(bool(records), "Cutout has no source-value bindings.")
        return records

    def coordinate(self, alias, rows, columns, artist, axis):
        self.bind(
            alias,
            rows,
            columns,
            artist.get_xdata if axis == "x" else artist.get_ydata,
            role=f"{axis}_coordinate",
        )

    def label(self, alias, row, column, artist, formatter=str):
        self.bind(
            alias, [row], [column], artist.get_text, lambda f: formatter(f.iloc[0, 0]), "label"
        )


def grid(ax, axis="y"):
    ax.set_axisbelow(True)
    ax.grid(axis=axis, color="#E2E6E9", linewidth=0.8)


def donor_labels(ax, trace, alias, rows):
    frame = trace.tables[alias].loc[rows]
    ax.set_yticks(range(len(frame)), labels=frame.donor.tolist())
    for row, artist in zip(rows, ax.get_yticklabels(), strict=True):
        trace.label(alias, row, "donor", artist)


def compartment_points(ax, x, y, comp):
    return ax.plot(
        x,
        y,
        linestyle="none",
        marker="o" if comp == "CD4" else "s",
        color=COLORS[comp],
        markerfacecolor=COLORS[comp] if comp == "CD4" else "white",
    )[0]


def key(fig, rect, entries):
    ax = fig.add_axes(rect)
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.set_axis_off()
    for x, label, marker, color, face in entries:
        ax.plot(x, 0.5, marker=marker, color=color, markerfacecolor=face, linestyle="none")
        ax.text(x + 0.03, 0.5, label, va="center", fontsize=24)


def design(fig, t):
    samples, depth = t.tables["rp_samples"], t.tables["rp_depth"]
    donors = sorted(samples.donor.unique())
    timeline = fig.add_axes((0.095, 0.20, 0.40, 0.67))
    coverage = fig.add_axes((0.655, 0.20, 0.30, 0.67))
    key(
        fig,
        (0.095, 0.90, 0.75, 0.09),
        [
            (0.025, "CD4", "o", COLORS["CD4"], COLORS["CD4"]),
            (0.23, "CD8", "s", COLORS["CD8"], "white"),
        ],
    )
    for y, donor in enumerate(donors):
        for comp, offset in [("CD4", -0.14), ("CD8", 0.14)]:
            rows = samples[samples.donor.eq(donor) & samples.compartment.eq(comp)].sort_values(
                "visit"
            )
            line = timeline.plot(
                rows.elapsed_years, np.full(len(rows), y + offset), color=COLORS[comp], alpha=0.5
            )[0]
            t.coordinate("rp_samples", rows.index, ["elapsed_years"], line, "x")
            marker = compartment_points(
                timeline, rows.elapsed_years, np.full(len(rows), y + offset), comp
            )
            t.coordinate("rp_samples", rows.index, ["elapsed_years"], marker, "x")
            row = depth[depth.donor.eq(donor) & depth.compartment.eq(comp)].iloc[0]
            bar = coverage.barh(
                y + offset,
                row.qualified_fraction_source,
                height=0.25,
                facecolor=COLORS[comp] if comp == "CD4" else "white",
                edgecolor=COLORS[comp],
                hatch=None if comp == "CD4" else "///",
            )[0]
            t.bind(
                "rp_depth",
                [row.name],
                ["qualified_fraction_source"],
                bar.get_width,
                lambda f: f.iloc[0, 0],
                "bar_width",
            )
    donor_rows = [samples[samples.donor.eq(donor)].index[0] for donor in donors]
    for ax in [timeline, coverage]:
        donor_labels(ax, t, "rp_samples", donor_rows)
        ax.set_ylim(5.6, -0.6)
        ax.spines["left"].set_visible(False)
        grid(ax, "x")
    timeline.set(xlim=(-0.6, 26), xticks=[0, 10, 20, 25], xlabel="Years from first visit")
    coverage.set(xlim=(0, 1), xticks=[0, 0.5, 1], xlabel="Qualified / source reads")
    coverage.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))


def cancellation(fig, t):
    source = t.tables["rp_controls"]
    key(
        fig,
        (0.085, 0.90, 0.85, 0.09),
        [
            (0.02, "CD4 RFU", "o", COLORS["CD4"], COLORS["CD4"]),
            (0.25, "CD8 RFU", "s", COLORS["CD8"], "white"),
            (0.49, "Matched maps: median / range", "D", GRAY, "white"),
        ],
    )
    for j, comp in enumerate(["CD4", "CD8"]):
        ax = fig.add_axes((0.095 + j * 0.50, 0.18, 0.37, 0.59))
        part = source[source.compartment.eq(comp)].sort_values("donor")
        for i, row in enumerate(part.itertuples()):
            connector = ax.plot(
                [row.observed_cancellation, row.control_median],
                [i - 0.15, i + 0.15],
                color="#C4CDD3",
            )[0]
            t.coordinate(
                "rp_controls",
                [row.Index],
                ["observed_cancellation", "control_median"],
                connector,
                "x",
            )
            interval = ax.plot(
                [row.control_q025, row.control_q975],
                [i + 0.15] * 2,
                color=COLORS[comp],
                linewidth=3,
            )[0]
            t.coordinate(
                "rp_controls", [row.Index], ["control_q025", "control_q975"], interval, "x"
            )
            median = ax.plot(
                row.control_median, i + 0.15, "D", color=COLORS[comp], markerfacecolor="white"
            )[0]
            t.coordinate("rp_controls", [row.Index], ["control_median"], median, "x")
            observed = compartment_points(ax, [row.observed_cancellation], [i - 0.15], comp)
            t.coordinate("rp_controls", [row.Index], ["observed_cancellation"], observed, "x")
        donor_labels(ax, t, "rp_controls", part.index.tolist())
        ax.set(
            xlim=(0, 0.7),
            xticks=[0, 0.2, 0.4, 0.6],
            ylim=(5.6, -0.6),
            xlabel="Aggregation cancellation",
        )
        title = ax.set_title(comp, color=COLORS[comp], fontweight="bold", pad=20)
        t.label("rp_controls", part.index[0], "compartment", title)
        ax.spines["left"].set_visible(False)
        grid(ax, "x")


def persistence(fig, t):
    frame = t.tables["rp_persistence"]
    metrics = ["persistent_fraction_union", "fraction_persistent_without_shared_clone"]
    labels = ["Persistent RFUs /\nunion detected RFUs", "No shared receptor /\npersistent RFUs"]
    for j, (metric, label) in enumerate(zip(metrics, labels, strict=True)):
        ax = fig.add_axes((0.13 + j * 0.51, 0.20, 0.31, 0.73))
        for i, donor in enumerate(sorted(frame.donor.unique())):
            part = frame[frame.donor.eq(donor)].sort_values("compartment")
            x = np.array([0, 1]) + (i - 2.5) * 0.035  # Display separation only.
            line = ax.plot(x, part[metric], color="#AAB5BD", alpha=0.65)[0]
            t.coordinate("rp_persistence", part.index, [metric], line, "y")
            for pos, row in zip(x, part.itertuples(), strict=True):
                marker = compartment_points(ax, [pos], [getattr(row, metric)], row.compartment)
                t.coordinate("rp_persistence", [row.Index], [metric], marker, "y")
        ax.set(
            xlim=(-0.4, 1.4),
            ylim=(0, 1.03),
            xticks=[0, 1],
            xticklabels=["CD4", "CD8"],
            yticks=[0, 0.5, 1],
            ylabel=label,
        )
        for comp, artist in zip(["CD4", "CD8"], ax.get_xticklabels(), strict=True):
            t.label(
                "rp_persistence", frame[frame.compartment.eq(comp)].index[0], "compartment", artist
            )
            artist.set_color(COLORS[comp])
        ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        grid(ax)


def gse190905(fig, t):
    frame = t.tables["gse190905_pairs"]
    placements = [
        ("all_T", (0.085, 0.19, 0.37, 0.67)),
        ("CD4", (0.61, 0.25, 0.145, 0.52)),
        ("CD8", (0.835, 0.25, 0.145, 0.52)),
    ]
    for subset, rect in placements:
        ax = fig.add_axes(rect)
        part = frame[frame.subset.eq(subset)].sort_values("donor")
        color = COLORS.get(subset, INK)
        for row in part.itertuples():
            line = ax.plot([0, 1], [row.d_clone, row.d_group], color=color, alpha=0.35)[0]
            t.coordinate("gse190905_pairs", [row.Index], ["d_clone", "d_group"], line, "y")
            marker = ax.plot(
                [0, 1],
                [row.d_clone, row.d_group],
                linestyle="none",
                marker="s" if subset == "CD8" else "o",
                markerfacecolor="white" if subset == "CD8" else color,
                color=color,
            )[0]
            t.coordinate("gse190905_pairs", [row.Index], ["d_clone", "d_group"], marker, "y")
        ax.set(
            xlim=(-0.4, 1.4),
            ylim=(0, 1.04),
            yticks=[0, 0.5, 1],
            xticks=[0, 1],
            xticklabels=["Receptor TV", "RFU TV"] if subset == "all_T" else ["Receptor", "RFU"],
        )
        title = ax.set_title(
            "Total T" if subset == "all_T" else subset,
            fontsize=32 if subset == "all_T" else 28,
            color=color,
            fontweight="bold",
            pad=18,
        )
        t.label(
            "gse190905_pairs",
            part.index[0],
            "subset",
            title,
            lambda value: "Total T" if value == "all_T" else value,
        )
        if subset == "all_T":
            ax.set_ylabel("Total variation (TV)")
        grid(ax)


def coverage_matrix(fig, t):
    source = t.tables["gse280982_visits"]
    groups = (
        source[["donor", "compartment"]].drop_duplicates().sort_values(["donor", "compartment"])
    )
    ax = fig.add_axes((0.23, 0.25, 0.61, 0.71))
    cmap = plt.get_cmap("Blues")
    for y, group in enumerate(groups.itertuples()):
        for visit in [1, 2, 3]:
            row = source[
                source.donor.eq(group.donor)
                & source.compartment.eq(group.compartment)
                & source.visit.eq(visit)
            ].iloc[0]
            available = row.source_available
            color = cmap(0.15 + 0.70 * row.coverage) if available else "#F1F2F3"
            cell = patches.Rectangle(
                (visit - 1.5, y - 0.5),
                1,
                1,
                facecolor=color,
                edgecolor="white",
                linewidth=3,
                hatch=None if available else "///",
            )
            ax.add_patch(cell)
            if available:
                t.bind(
                    "gse280982_visits",
                    [row.name],
                    ["coverage"],
                    cell.get_facecolor,
                    lambda f: cmap(0.15 + 0.70 * f.iloc[0, 0]),
                    "coverage_color",
                )
            t.bind(
                "gse280982_visits",
                [row.name],
                ["source_available"],
                cell.get_hatch,
                lambda f: None if f.iloc[0, 0] else "///",
                "missing_hatch",
            )

            def formatter(f):
                return (
                    f"{int(f.qualified_cells.iloc[0]):,}/{int(f.primary_TRB_cells.iloc[0]):,}"
                    if f.source_available.iloc[0]
                    else "Missing"
                )

            label = formatter(source.loc[[row.name]])
            artist = ax.text(
                visit - 1,
                y,
                label,
                va="center",
                ha="center",
                fontsize=24,
                color="white" if available and row.coverage > 0.65 else INK,
            )
            t.bind(
                "gse280982_visits",
                [row.name],
                ["source_available", "qualified_cells", "primary_TRB_cells"],
                artist.get_text,
                formatter,
                "visit_counts_or_missing",
            )
    ax.set(
        xlim=(-0.5, 2.5),
        ylim=(5.5, -0.5),
        xticks=[0, 1, 2],
        xticklabels=["Pre", "Last radiation\nday", "6 weeks"],
        yticks=range(6),
        yticklabels=[f"{r.donor} {r.compartment}" for r in groups.itertuples()],
        xlabel="Qualified / primary-TRB cells",
    )
    for row, artist in zip(groups.index, ax.get_yticklabels(), strict=True):
        t.bind(
            "gse280982_visits",
            [row],
            ["donor", "compartment"],
            artist.get_text,
            lambda f: f"{f.donor.iloc[0]} {f.compartment.iloc[0]}",
            "visit_row_label",
        )
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="x", length=0)
    visit_labels = {1: "Pre", 2: "Last radiation\nday", 3: "6 weeks"}
    for visit, artist in zip([1, 2, 3], ax.get_xticklabels(), strict=True):
        t.label(
            "gse280982_visits",
            source[source.visit.eq(visit)].index[0],
            "visit",
            artist,
            lambda value: visit_labels[value],
        )
    cb = fig.colorbar(
        plt.cm.ScalarMappable(
            norm=Normalize(0, 1), cmap=ListedColormap(cmap(np.linspace(0.15, 0.85, 256)))
        ),
        cax=fig.add_axes((0.88, 0.28, 0.02, 0.65)),
        ticks=[0, 0.5, 1],
        format=PercentFormatter(1, decimals=0),
    )
    cb.outline.set_visible(False)
    cb.solids.set_rasterized(False)


INTERVALS = {
    "pre_to_last_radiation": (0, "pre–radiation"),
    "last_radiation_to_6weeks": (1, "radiation–6 weeks"),
    "pre_to_6weeks": (2, "pre–6 weeks"),
}


def interval_tv(fig, t):
    frame = t.tables["gse280982_pairs"]
    supported = frame[frame.status.eq("analyzed")]
    rows = sorted(
        supported.index,
        key=lambda i: (
            frame.loc[i, "compartment"],
            frame.loc[i, "donor"],
            INTERVALS[frame.loc[i, "interval"]][0],
        ),
    )
    ax = fig.add_axes((0.42, 0.15, 0.54, 0.69))
    key(
        fig,
        (0.42, 0.89, 0.54, 0.09),
        [(0.03, "Receptor TV", "o", GRAY, "white"), (0.58, "RFU TV", "D", RFU, RFU)],
    )
    for y, idx in enumerate(rows):
        row = frame.loc[idx]
        line = ax.plot([row.d_clone, row.d_group], [y, y], color="#A8B3BB", linewidth=2.3)[0]
        t.coordinate("gse280982_pairs", [idx], ["d_clone", "d_group"], line, "x")
        for col, marker, color, face in [
            ("d_clone", "o", GRAY, "white"),
            ("d_group", "D", RFU, RFU),
        ]:
            point = ax.plot(row[col], y, marker=marker, color=color, markerfacecolor=face)[0]
            t.coordinate("gse280982_pairs", [idx], [col], point, "x")

    def formatter(f):
        return f"{f.donor.iloc[0]} {f.compartment.iloc[0]} · {INTERVALS[f.interval.iloc[0]][1]}"

    ax.set(
        xlim=(0, 1.025),
        xticks=[0, 0.5, 1],
        ylim=(7.5, -0.5),
        yticks=range(8),
        yticklabels=[formatter(frame.loc[[i]]) for i in rows],
        xlabel="Total variation (TV)",
    )
    for idx, artist in zip(rows, ax.get_yticklabels(), strict=True):
        t.bind(
            "gse280982_pairs",
            [idx],
            ["donor", "compartment", "interval"],
            artist.get_text,
            formatter,
            "interval_label",
        )
    grid(ax, "x")


def regulatory_counts(fig, t):
    frame = t.tables["regulatory_counts"]
    categories = ["published_RFU_QTL", "conditional_eQTL", "conditional_caQTL", "same_variant_both"]
    rows = [frame[frame.category.eq(category)].index[0] for category in categories]
    for idx, x, size in [
        (rows[0], 0.13, 64),
        (rows[1], 0.425, 56),
        (rows[2], 0.555, 56),
        (rows[3], 0.86, 64),
    ]:
        text = fig.text(
            x,
            0.68,
            f"{int(frame.loc[idx, 'variants']):,}",
            ha="center",
            va="center",
            fontsize=size,
            fontweight="bold",
        )
        t.label("regulatory_counts", idx, "variants", text, lambda value: f"{int(value):,}")
    for x in [0.28, 0.71]:
        fig.text(x, 0.68, "→", ha="center", va="center", fontsize=40, color=GRAY)
    fig.text(0.49, 0.68, "/", ha="center", va="center", fontsize=36, color=GRAY)
    for x, label in [
        (0.13, "RFU-QTL\nvariants"),
        (0.49, "Conditional eQTL /\ncaQTL overlaps"),
        (0.86, "Both layers\n(same variant)"),
    ]:
        fig.text(x, 0.34, label, ha="center", va="top", fontsize=24)


def regulatory_matrix(fig, t):
    frame = t.tables["regulatory_matrix"]
    ax = fig.add_axes((0.025, 0.04, 0.95, 0.73))
    ax.set(xlim=(0, 1), ylim=(5.6, -0.6))
    ax.set_axis_off()
    cols = [
        (0.005, "Variant (chr:pos:ref:alt)", "variant"),
        (0.335, "eQTL target", "eqtl_target"),
        (0.615, "Target class", "target_class"),
        (0.815, "caQTL\nrecords", "caqtl_records"),
        (0.96, "RFUs", "rfu_count"),
    ]
    for x, label, _ in cols:
        ax.text(
            x,
            -0.88,
            label,
            ha="center" if x > 0.8 else "left",
            va="bottom",
            fontsize=24,
            fontweight="bold",
        )
    for y, row in enumerate(frame.itertuples()):
        ax.add_patch(
            patches.Rectangle(
                (0, y - 0.43),
                1,
                0.86,
                facecolor="#F1F4F6" if y % 2 == 0 else "white",
                edgecolor="none",
            )
        )
        for x, _, col in cols:

            def formatter(value, column=col):
                if column == "target_class":
                    return {"direct_TCR_gene": "TCR", "non_TCR_gene": "Non-TCR"}[value]
                return str(int(value)) if column in {"caqtl_records", "rfu_count"} else str(value)

            artist = ax.text(
                x,
                y,
                formatter(getattr(row, col)),
                ha="center" if x > 0.8 else "left",
                va="center",
                fontsize=24,
                fontweight="bold" if col == "target_class" else "normal",
            )
            t.label("regulatory_matrix", row.Index, col, artist, formatter)


DRAW = {
    "rp1_14_design_coverage": design,
    "rp1_14_cancellation_controls": cancellation,
    "rp1_14_persistence_turnover": persistence,
    "gse190905_paired_tv": gse190905,
    "gse280982_visit_coverage": coverage_matrix,
    "gse280982_interval_tv": interval_tv,
    "regulatory_overlap_counts": regulatory_counts,
    "regulatory_six_variant_matrix": regulatory_matrix,
}


def read_sources(workspace, spec):
    roots, upstream, inputs, tables = {}, {}, {}, {}
    for name, record in spec["upstream"].items():
        root = workspace / record["directory"]
        completion = root / "completion.json"
        require(
            sha256(completion) == record["completion_sha256"],
            f"Changed completed figure pack: {root}",
        )
        source = json.loads(completion.read_text())
        require(source["status"] == "complete", f"Incomplete figure pack: {root}")
        for relative, checksum in source["outputs"].items():
            require(
                sha256(safe_child(root, relative)) == checksum,
                f"Changed frozen output: {root / relative}",
            )
        roots[name] = root
        upstream[name] = {
            "completion": str(completion),
            "sha256": sha256(completion),
            "verified_outputs": len(source["outputs"]),
        }
    for alias, (stage, relative) in spec["sources"].items():
        path = safe_child(roots[stage], relative)
        inputs[alias] = {
            "path": str(path),
            "sha256": sha256(path),
            "pack": stage,
            "copy": f"sources/{alias}.tsv",
        }
        tables[alias] = pd.read_csv(path, sep="\t", float_precision="round_trip")
    return tables, inputs, upstream


def export(fig, out, name, spec, asset, trace):
    audit = layout_audit(fig, spec["minimum_font_pt"])
    require(fig._suptitle is None, "Cutout acquired an overall title.")
    values = trace.check()
    json_save(out / "manifests" / f"{name}_value_checks.json", values)
    files = {}
    for extension in ["pdf", "svg", "png"]:
        path = out / "assets" / f"{name}.{extension}"
        metadata = {"Creator": "scRFU frozen poster-asset extraction"}
        if extension == "pdf":
            metadata.update(CreationDate=None, ModDate=None)
        elif extension == "svg":
            metadata["Date"] = None
        fig.savefig(path, dpi=spec["png_dpi"], metadata=metadata)
        files[str(path.relative_to(out))] = sha256(path)
    fig.savefig(out / "preview" / f"{name}.png", dpi=110)
    plt.close(fig)
    (out / "captions" / f"{name}.md").write_text(asset["caption"] + "\n")
    return {
        "name": name,
        "source_tables": asset["sources"],
        "files_sha256": files,
        "caption_file": f"captions/{name}.md",
        "size_inches": asset["size_inches"],
        "placement_mm": [v * 25.4 for v in asset["size_inches"]],
        "png_dpi": spec["png_dpi"],
        "verified_artist_bindings": len(values),
        "value_checks": f"manifests/{name}_value_checks.json",
        **audit,
    }


def contact_sheet(out, spec):
    sheet = Image.new("RGB", (2080, 2360), "#E3E8EC")
    text = ImageDraw.Draw(sheet)
    font = ImageFont.truetype(font_manager.findfont("DejaVu Sans"), 30)
    for i, (name, asset) in enumerate(spec["assets"].items()):
        left, top = 20 + i % 2 * 1040, 20 + i // 2 * 590
        text.text((left, top), asset["contact_label"], fill=INK, font=font)
        with Image.open(out / "preview" / f"{name}.png") as image:
            thumbnail = image.convert("RGB")
            thumbnail.thumbnail((1000, 500))
            sheet.paste(thumbnail, (left, top + 55))
    sheet.save(out / "contact_sheet.png")


def validate(out, inputs, completion):
    for relative, checksum in completion["outputs"].items():
        require(sha256(safe_child(out, relative)) == checksum, f"Changed asset: {relative}")
    for record in inputs.values():
        require(
            sha256(Path(record["path"])) == record["sha256"] == sha256(out / record["copy"]),
            "Frozen TSV copy differs.",
        )
    manifest = json.loads((out / "manifests/asset_manifest.json").read_text())
    probes = 0
    for asset in manifest["exports"]:
        require(
            asset["minimum_font_pt"] >= 24
            and not asset["text_outside_canvas"]
            and not asset["text_overlaps"],
            "Asset typography failed.",
        )
        require((out / asset["caption_file"]).is_file(), "Missing caption.")
        checks = json.loads((out / asset["value_checks"]).read_text())
        probes += len(checks)
        # Re-read unchanged TSV copies with round-trip precision, independent of artists.
        for check in checks:
            source = pd.read_csv(
                out / inputs[check["table"]]["copy"], sep="\t", float_precision="round_trip"
            )
            values = source.loc[
                [line - 2 for line in check["tsv_line_numbers"]], check["columns"]
            ].to_numpy()
            require(
                plain(values) == check["frozen_values"],
                "Recorded artist reference differs from frozen table.",
            )
        for relative, checksum in asset["files_sha256"].items():
            path = out / relative
            require(sha256(path) == checksum, "Asset manifest output hash differs.")
            if path.suffix == ".svg":
                require(
                    not any(e.tag.endswith("}image") for e in ET.parse(path).getroot().iter()),
                    "Rasterized SVG content.",
                )
            elif path.suffix == ".pdf":
                require(
                    len(
                        subprocess.check_output(["pdfimages", "-list", str(path)], text=True)
                        .strip()
                        .splitlines()
                    )
                    == 2,
                    "Rasterized PDF content.",
                )
                info = subprocess.check_output(["pdfinfo", str(path)], text=True)
                require(
                    next(line for line in info.splitlines() if line.startswith("Pages:")).split()[
                        -1
                    ]
                    == "1",
                    "Multi-page cutout.",
                )
            else:
                with Image.open(path) as image:
                    require(
                        image.size == tuple(round(v * 300) for v in asset["size_inches"]),
                        "PNG dimensions differ.",
                    )
                    require(
                        all(abs(v - 300) < 0.1 for v in image.info.get("dpi", (0, 0))),
                        "PNG resolution differs.",
                    )
                    image.verify()
    require(len(manifest["exports"]) == 8, "Incomplete cutout pack.")
    return {
        "status": "passed",
        "cutouts": 8,
        "asset_files": 24,
        "source_copies": len(inputs),
        "verified_artist_bindings": probes,
        "all_values_match_frozen_tables": True,
        "minimum_font_pt": 24,
        "vector_pdfs_and_svgs": True,
        "biological_analysis_or_new_results": False,
    }


def run(args):
    # Completed proof is metadata about one export, not an input to rendering.
    spec = {
        key: value for key, value in json.loads(MANIFEST.read_text()).items() if key != "completed"
    }
    workspace = (args.workspace or Path(spec["workspace"])).resolve()
    out = (args.out or workspace / spec["output_directory"]).resolve()
    require(not out.is_relative_to(REPO), "Generated assets/source copies must remain outside Git.")
    for source in spec["upstream"].values():
        require(
            not out.is_relative_to(workspace / source["directory"]),
            "Cannot write inside a completed figure pack.",
        )
    tables, inputs, upstream = read_sources(workspace, spec)
    identity = {
        "spec": spec,
        "inputs": inputs,
        "upstream": upstream,
        "python": sys.version,
        "executable": sys.executable,
        "versions": {
            name: importlib.metadata.version(name)
            for name in ["matplotlib", "pandas", "numpy", "Pillow"]
        },
        "code_sha256": {
            str(path.relative_to(REPO)): sha256(path)
            for path in [
                Path(__file__),
                Path(__file__).with_name("rp1_14_reusable_figures.py"),
                Path(__file__).with_name("_rp1_14_figure_drawing.py"),
            ]
        },
    }
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    marker = out / "completion.json"
    if marker.exists():
        completion = json.loads(marker.read_text())
        require(
            completion["status"] == "complete" and completion["fingerprint"] == fingerprint,
            "Completed extraction identity changed; use a new versioned directory.",
        )
        print(json.dumps(validate(out, inputs, completion), sort_keys=True))
        print("Verified without rewriting assets or source figures.")
        return
    require(not args.verify_only, "No completed asset pack to verify.")
    require(not out.exists() or not any(out.iterdir()), "Output nonempty; use a new directory.")
    for name in ["assets", "sources", "captions", "manifests", "preview"]:
        (out / name).mkdir(parents=True, exist_ok=True)
    for record in inputs.values():
        shutil.copyfile(record["path"], out / record["copy"])
    command = shlex.join(
        [
            sys.executable,
            "-m",
            "manuscript.scripts.apbc2026_poster_assets",
            "--workspace",
            str(workspace),
            "--out",
            str(out),
        ]
    )
    exports = []
    with plt.rc_context(STYLE):
        for name, asset in spec["assets"].items():
            trace = Trace(tables)
            fig = plt.figure(figsize=asset["size_inches"])
            DRAW[name](fig, trace)
            exports.append(export(fig, out, name, spec, asset, trace))
    contact_sheet(out, spec)
    json_save(
        out / "manifests/asset_manifest.json",
        {
            **identity,
            "fingerprint": fingerprint,
            "command": command,
            "exports": exports,
            "new_results": [],
            "source_tables_copied_byte_for_byte": True,
        },
    )
    (out / "captions.md").write_text(
        "# Poster cutout captions\n\n"
        + "\n\n".join(f"**{name}**: {asset['caption']}" for name, asset in spec["assets"].items())
        + "\n"
    )
    (out / "README.md").write_text(
        "# APBC 2026 poster assets v1\n\nEight standalone cutouts from the two completed reusable figure packs. No analysis, new summary, ratio, median, resampling or conceptual cartoon was produced.\n\n`assets/`: vector PDF/SVG and 300-dpi PNG; `contact_sheet.png`: visual selection sheet; `captions/`: one-line captions; `sources/`: unchanged TSV copies; `manifests/`: source/output hashes and checks of actual plotted coordinates, labels and colors. No overall title or panel letter appears on a cutout.\n\nAll cutouts use a white background and a minimum 24-point font at their native size. The asset manifest records size in inches and millimeters for placement in PowerPoint/BioRender. Keep native size or enlarge; shrinking also shrinks labels. The contact sheet is a selection preview, not an A0 print assembly.\n\nReproduce from the scRFU repository root:\n\n```bash\n"
        + command
        + "\n```\n\nAppend `--verify-only` to validate without exporting. Normal repeats also verify without rewriting. Changed inputs, code or software require a new directory. All 252 outputs in the two completed figure packs are hash-verified, and all nine source TSVs are copied byte-for-byte. Manuscript figures remain untouched.\n\nEvery measured numerical coordinate, control range, matrix count, coverage color and data label is linked to its original TSV line and column and checked after rendering. Percentage ticks and color scales are display formatting only. No new donor summary bars are computed; the matched-map medians and ranges are already present in the frozen table. Regulatory arrows indicate set-overlap counts, not a biological mechanism or formal colocalization.\n"
    )
    initial = {
        "outputs": {
            str(path.relative_to(out)): sha256(path)
            for path in sorted(out.rglob("*"))
            if path.is_file()
        }
    }
    report = validate(out, inputs, initial)
    report["upstream_outputs_verified"] = sum(
        record["verified_outputs"] for record in upstream.values()
    )
    report["completed_figure_packs_unchanged"] = True
    json_save(out / "manifests/validation.json", report)
    outputs = {
        str(path.relative_to(out)): sha256(path)
        for path in sorted(out.rglob("*"))
        if path.is_file()
    }
    json_save(marker, {"status": "complete", "fingerprint": fingerprint, "outputs": outputs})
    print(json.dumps(report, sort_keys=True))
    print(f"Poster cutouts complete: {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--verify-only", action="store_true")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
