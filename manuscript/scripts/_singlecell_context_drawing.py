"""Publication-style context plots from frozen figure tables (no statistics)."""

from __future__ import annotations

import hashlib
import json
import subprocess
import xml.etree.ElementTree as ET

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager, patches, transforms
from matplotlib.textpath import TextPath
from PIL import Image, ImageDraw, ImageFont

from ._rp1_14_figure_drawing import GRAY, INK, layout_audit, style
from .rp1_14_reusable_figures import json_save, require, sha256

PURPLE = "#72549A"
SOURCE_LABELS = {
    "naive thymus-derived CD4-positive, alpha-beta T cell": ("CD4 naive", "#86BDDC"),
    "central memory CD4-positive, alpha-beta T cell": ("CD4 central memory", "#0072B2"),
    "CD4-positive, alpha-beta memory T cell": ("CD4 memory", "#224F73"),
    "effector memory CD4-positive, alpha-beta T cell": ("CD4 effector memory", "#4B8DAB"),
    "effector memory CD4-positive, alpha-beta T cell, terminally differentiated": (
        "CD4 terminal effector memory",
        "#102D43",
    ),
    "CD4-positive, CD25-positive, alpha-beta regulatory T cell": ("CD4 regulatory", "#83D0C0"),
    "naive thymus-derived CD8-positive, alpha-beta T cell": ("CD8 naive", "#E9BD82"),
    "central memory CD8-positive, alpha-beta T cell": ("CD8 central memory", "#CC913F"),
    "CD8-positive, alpha-beta memory T cell": ("CD8 memory", "#C15A00"),
    "effector memory CD8-positive, alpha-beta T cell": ("CD8 effector memory", "#E78953"),
    "effector memory CD8-positive, alpha-beta T cell, terminally differentiated": (
        "CD8 terminal effector memory",
        "#7D3D0F",
    ),
    "gamma-delta T cell": ("Gamma-delta T", "#A16EA7"),
    "mucosal invariant T cell": ("Mucosal invariant T", "#705696"),
    "germinal center T cell": ("Germinal center T", "#458258"),
    "group 3 innate lymphoid cell": ("Group 3 innate lymphoid", "#727D49"),
    "innate lymphoid cell": ("Innate lymphoid", "#8B8F98"),
    "CD16-positive, CD56-dim natural killer cell, human": (
        "CD16+ / CD56-dim natural killer",
        "#4B4C51",
    ),
}


class Checks:
    def __init__(self):
        self.records = []

    def array(self, table, role, expected, actual):
        expected, actual = (
            np.asarray(expected, dtype=np.float64),
            np.asarray(actual, dtype=np.float64),
        )
        require(np.array_equal(expected, actual), f"Plotted value mismatch: {table} {role}")
        self.records.append(
            {
                "table": table,
                "role": role,
                "shape": list(expected.shape),
                "sha256_float64": hashlib.sha256(expected.tobytes()).hexdigest(),
                "exact_match": True,
            }
        )

    def text(self, table, role, expected, artist):
        require(artist.get_text() == expected, f"Plotted text mismatch: {role}")
        self.records.append({"table": table, "role": role, "text": expected, "exact_match": True})


def frame(fig, rect, title=None, subtitle=None):
    x, y, w, h = rect
    if title:
        fig.text(x, y + h, title, ha="left", va="top", fontweight="bold")
    if subtitle:
        fig.text(x, y + h - 0.045, subtitle, ha="left", va="top", color=GRAY)


def umap(fig, rect, tables, checks, view, summary, *, title=True):
    points = tables["A_umap_cells"]
    x, y, w, h = rect
    ax = fig.add_axes((x + w * 0.055, y + h * 0.10, w * 0.93, h * 0.72))
    coordinates = points[["umap_1", "umap_2"]].to_numpy()
    limits = [(coordinates[:, i].min() - 0.6, coordinates[:, i].max() + 0.6) for i in range(2)]
    size = 2.2 * (fig.get_figwidth() / 8) ** 1.3

    def cloud(part, color, role, multiplier=1):
        artist = ax.scatter(
            part.umap_1,
            part.umap_2,
            s=size * multiplier,
            c=color,
            alpha=0.75 if view == "annotation" else 0.85,
            linewidths=0,
            rasterized=True,
        )
        checks.array(
            "A_umap_cells", f"{view}: {role}", part[["umap_1", "umap_2"]], artist.get_offsets()
        )

    if view == "annotation":
        require(set(points.cell_type) <= set(SOURCE_LABELS), "Unmapped source annotation")
        for label, (_, color) in SOURCE_LABELS.items():
            cloud(points[points.cell_type.eq(label)], color, label)
        heading, sub = "Source cell types", "Original atlas annotations"
    else:
        if view == "assignment":
            active = points.assignment_status.eq("qualified")
            heading, sub = (
                "Qualified RFU assignment",
                f"{summary['qualified_cells']:,} / {summary['atlas_cells']:,} atlas cells",
            )
        else:
            active = points.rfu.eq(view)
            support = tables["B_rfu_support"].set_index("rfu").loc[view]
            heading, sub = view, f"{support.cells:,} cells · {support.donors} donors"
        cloud(points[~active], "#D9DDE1", "context")
        cloud(
            points[active],
            PURPLE,
            "qualified" if view == "assignment" else view,
            1 if view == "assignment" else 4,
        )
    ax.set(
        xlim=limits[0],
        ylim=limits[1],
        xticks=[],
        yticks=[],
        xlabel="UMAP 1",
        ylabel="UMAP 2",
        aspect="equal",
    )
    for spine in ax.spines.values():
        spine.set_visible(False)
    if title:
        frame(fig, rect, heading, sub)
    return ax


def annotation_key(fig, rect, tables):
    x, y, w, h = rect
    fig.text(x, y + h, "Source annotation key", va="top", fontweight="bold")
    for i, (_source, (label, color)) in enumerate(SOURCE_LABELS.items()):
        ypos = y + h * (0.9 - i * 0.053)
        fig.patches.append(
            patches.Circle(
                (x + w * 0.015, ypos), radius=0.0045, transform=fig.transFigure, color=color
            )
        )
        fig.text(x + w * 0.045, ypos, label, va="center")


def aa_color(aa):
    for letters, color in [
        ("KRH", "#C15A00"),
        ("DE", "#0072B2"),
        ("STNQ", "#258570"),
        ("GP", PURPLE),
        ("C", "#9A791B"),
    ]:
        if aa in letters:
            return color
    return INK


def logos(fig, rect, tables, checks, labels, *, compact=False):
    x, y, w, h = rect
    support = tables["B_rfu_support"].set_index("rfu")
    frequencies = tables["B_logo_frequencies"]
    font = font_manager.FontProperties(family="DejaVu Sans", weight="bold")
    for i, label in enumerate(labels):
        row = support.loc[label]
        ry = y + h - (i + 1) * h / len(labels)
        rh = h / len(labels)
        heading = (
            f"{label}   ·   {row.logo_unique_aa} unique sequences at length {row.modal_length}"
        )
        artist = fig.text(x + w * 0.05, ry + rh * 0.96, heading, fontweight="bold", va="top")
        checks.text("B_rfu_support", f"{label}: modal subset", heading, artist)
        sub = f"Family support: {row.cells:,} cells · {row.unique_aa} unique amino-acid sequences · {row.donors} donors"
        artist = fig.text(x + w * 0.05, ry + rh * 0.79, sub, color=GRAY, va="top")
        checks.text("B_rfu_support", f"{label}: support", sub, artist)
        ax = fig.add_axes((x + w * 0.06, ry + rh * 0.20, w * 0.92, rh * 0.43))
        frame_values = frequencies[frequencies.rfu.eq(label)]
        glyph_values = []
        for position, part in frame_values.groupby("position", sort=True):
            base = 0
            for r in part.sort_values(["frequency", "amino_acid"]).itertuples():
                if r.frequency == 0:
                    continue
                path = TextPath((0, 0), r.amino_acid, size=1, prop=font)
                bounds = path.get_extents()
                transform = (
                    transforms.Affine2D()
                    .translate(-bounds.x0, -bounds.y0)
                    .scale(0.88 / bounds.width, r.frequency / bounds.height)
                    .translate(position - 0.44, base)
                )
                patch = patches.PathPatch(
                    path,
                    color=aa_color(r.amino_acid),
                    linewidth=0,
                    transform=transform + ax.transData,
                )
                ax.add_patch(patch)
                height = transform.transform_path(path).get_extents().height
                require(
                    abs(height - r.frequency) < 1e-14, "Logo glyph height differs from frequency"
                )
                glyph_values.append((position, ord(r.amino_acid), r.frequency))
                base += r.frequency
            require(abs(base - 1) < 1e-14, "Logo position does not sum to one")
        expected = frame_values[frame_values.frequency.gt(0)].sort_values(
            ["position", "frequency", "amino_acid"]
        )
        checks.array(
            "B_logo_frequencies",
            label,
            np.column_stack([expected.position, expected.amino_acid.map(ord), expected.frequency]),
            glyph_values,
        )
        ax.set(
            xlim=(0.5, row.modal_length + 0.5),
            ylim=(0, 1.02),
            yticks=[0, 1],
            xticks=sorted(set([1, 5, 10, int(row.modal_length)])),
            ylabel="Frequency",
        )
        if i == len(labels) - 1:
            ax.set_xlabel("CDR3 amino-acid position")
        ax.spines["left"].set_visible(False)
    if not compact:
        fig.text(
            x + w * 0.05,
            y - 0.065,
            "Sequence-family summaries · modal length only · no antigen inference",
            color=GRAY,
        )


def volcano(fig, rect, tables, checks, summary):
    x, y, w, h = rect
    frame(
        fig,
        rect,
        "Exploratory expression contrast",
        f"RFU3526 vs matched other RFUs · {summary['de_donors']} donors",
    )
    ax = fig.add_axes((x + w * 0.14, y + h * 0.22, w * 0.78, h * 0.61))
    data = tables["C_differential_expression"]
    data = data[data.test_status.eq("tested")].copy()
    # Show raw P spread while retaining the adjusted-significance color/count.
    data["volcano_y"] = -np.log10(data.p_value.astype(float))
    for significant, color in [(False, "#A4AEB7"), (True, PURPLE)]:
        part = data[data.q_lt_0_05.eq(significant)]
        coords = part[["mean_paired_log2cpm_difference", "volcano_y"]].to_numpy(dtype=float)
        artist = ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=color,
            s=8 * (fig.get_figwidth() / 8),
            alpha=0.7,
            linewidths=0,
            rasterized=True,
        )
        checks.array(
            "C_differential_expression",
            f"x=paired effect, y=-log10(raw P), q<0.05={significant}",
            coords,
            artist.get_offsets(),
        )
    extent = max(1, float(data.mean_paired_log2cpm_difference.abs().max()) * 1.12)
    ymax = max(2.3, float(data.volcano_y.max()) * 1.5)
    ax.axvline(0, color="#CAD1D6", linewidth=1, zorder=0)
    selected = data[data.label_gene].sort_values("mean_paired_log2cpm_difference")
    locations = [(-0.70, 0.88), (-0.25, 0.70), (0.25, 0.88), (0.70, 0.70)]
    for (_, row), (tx, ty) in zip(selected.iterrows(), locations, strict=True):
        ax.plot(
            [float(row.mean_paired_log2cpm_difference), extent * tx],
            [float(row.volcano_y), ymax * (ty - 0.025)],
            color=GRAY,
            linewidth=0.7,
        )
        artist = ax.text(extent * tx, ymax * ty, row.gene, ha="center", va="center")
        checks.text("C_differential_expression", f"top-gene label: {row.gene_id}", row.gene, artist)
    ax.set(
        xlim=(-extent, extent),
        ylim=(0, ymax),
        ylabel="−log₁₀ raw P",
        xlabel="Mean paired Δ log₂(CPM + 1)\nRFU3526 − matched background",
    )
    ax.locator_params(axis="x", nbins=5)
    ax.locator_params(axis="y", nbins=4)
    ax.tick_params(axis="x", pad=7)
    ax.tick_params(axis="y", pad=7)
    note = f"{summary['de_cells_per_group']} cells/group · {summary['de_cells_per_donor_min']}–{summary['de_cells_per_donor_max']} per donor/group"
    artist = fig.text(x + w * 0.14, y + h * 0.03, note, color=GRAY)
    checks.text("C_pseudobulk_samples", "contrast sample sizes", note, artist)
    note = f"{summary['de_genes_q_lt_0_05']:,} / {summary['de_genes_tested']:,} genes with q < 0.05"
    artist = ax.text(0.98, 0.97, note, ha="right", va="top", transform=ax.transAxes, color=GRAY)
    checks.text("C_differential_expression", "adjusted significance summary", note, artist)


def draw(name, mode, tables, summary, spec):
    # Dimensions are native placement sizes, not bounding-box-cropped canvases.
    sizes = {
        "A_umap_context": (24, 15),
        "B_rfu_logos": (18, 12),
        "C_exploratory_volcano": (13, 11),
        "D_context_block": (26, 18),
    }
    manuscript = {
        "A_umap_context": (9.6, 6),
        "B_rfu_logos": (7.2, 4.8),
        "C_exploratory_volcano": (5.2, 4.4),
        "D_context_block": (10.4, 7.2),
    }
    size = sizes[name] if mode == "poster" else manuscript[name]
    fig = plt.figure(figsize=size)
    checks = Checks()
    labels = spec["selected_rfus"]
    if name == "A_umap_context":
        views = ["annotation", "assignment", *labels]
        for i, view in enumerate(views):
            umap(
                fig,
                (0.025 + i % 3 * 0.33, 0.55 if i < 3 else 0.075, 0.30, 0.40),
                tables,
                checks,
                view,
                summary,
            )
        annotation_key(fig, (0.69, 0.075, 0.29, 0.40), tables)
        fig.text(
            0.025,
            0.012,
            f"Frozen Wells UMAP · {summary['displayed_umap_cells']:,} cells displayed · includes all cells from the three selected RFUs",
            color=GRAY,
        )
    elif name == "B_rfu_logos":
        logos(fig, (0.02, 0.11, 0.95, 0.87), tables, checks, labels)
    elif name == "C_exploratory_volcano":
        volcano(fig, (0.015, 0.02, 0.975, 0.95), tables, checks, summary)
    else:
        umap(fig, (0.025, 0.57, 0.40, 0.40), tables, checks, labels[0], summary)
        volcano(fig, (0.025, 0.01, 0.42, 0.52), tables, checks, summary)
        logos(fig, (0.48, 0.23, 0.50, 0.74), tables, checks, labels, compact=True)
        fig.text(0.515, 0.19, "Sequence families; no antigen inference", color=GRAY)
        # Lower-right region intentionally reserved for later BioRender composition.
    return fig, checks


def export(fig, out, mode, name, checks, minimum):
    audit = layout_audit(fig, minimum)
    files = {}
    for extension in ["pdf", "svg", "png"]:
        path = out / "figures" / mode / f"{name}.{extension}"
        metadata = {"Creator": "scRFU single-cell context figures"}
        if extension == "pdf":
            metadata.update(CreationDate=None, ModDate=None)
        elif extension == "svg":
            metadata["Date"] = None
        fig.savefig(path, dpi=300, metadata=metadata)
        files[str(path.relative_to(out))] = sha256(path)
    json_save(out / "manifests" / f"{mode}_{name}_value_checks.json", checks.records)
    preview = out / "preview" / f"{mode}_{name}.png"
    fig.savefig(preview, dpi=80 if mode == "poster" else 150)
    size = list(fig.get_size_inches())
    plt.close(fig)
    return {
        "figure": name,
        "mode": mode,
        "caption_file": "manifests/captions.md",
        "size_inches": size,
        "placement_mm": [v * 25.4 for v in size],
        "files_sha256": files,
        "source_tables": sorted({c["table"] for c in checks.records}),
        "value_checks": f"manifests/{mode}_{name}_value_checks.json",
        "artist_checks": len(checks.records),
        "dense_points_rasterized_dpi": 300,
        "text_axes_and_logos_vector": True,
        **audit,
    }


def render_all(out, spec, summary):
    tables = {
        p.stem: pd.read_csv(p, sep="\t", float_precision="round_trip", keep_default_na=False)
        for p in (out / "tables").glob("*.tsv")
    }
    annotation = pd.DataFrame(
        [
            {"source_cell_type": k, "display_label": v[0], "color": v[1]}
            for k, v in SOURCE_LABELS.items()
        ]
    )
    annotation.to_csv(out / "tables/A_annotation_key.tsv", sep="\t", index=False)
    exports = []
    for mode in ["manuscript", "poster"]:
        (out / "figures" / mode).mkdir(parents=True)
        (out / "preview").mkdir(exist_ok=True)
        settings = style(mode)
        minimum = spec[f"{mode}_min_font_pt"]
        settings.update(
            {
                key: minimum
                for key in [
                    "font.size",
                    "axes.labelsize",
                    "axes.titlesize",
                    "xtick.labelsize",
                    "ytick.labelsize",
                    "legend.fontsize",
                ]
            }
        )
        settings["svg.hashsalt"] = "singlecell-context-v1"
        with plt.rc_context(settings):
            for name in [
                "A_umap_context",
                "B_rfu_logos",
                "C_exploratory_volcano",
                "D_context_block",
            ]:
                print(f"Rendering {mode}/{name}", flush=True)
                fig, checks = draw(name, mode, tables, summary, spec)
                exports.append(export(fig, out, mode, name, checks, minimum))
    sheet = Image.new("RGB", (2100, 1700), "#E5E9EC")
    pen = ImageDraw.Draw(sheet)
    font = ImageFont.truetype(font_manager.findfont("DejaVu Sans"), 30)
    for i, name in enumerate(
        ["A_umap_context", "B_rfu_logos", "C_exploratory_volcano", "D_context_block"]
    ):
        left, top = 25 + i % 2 * 1050, 25 + i // 2 * 850
        pen.text((left, top), name.replace("_", " "), font=font, fill=INK)
        with Image.open(out / "preview" / f"poster_{name}.png") as source:
            thumb = source.convert("RGB")
            thumb.thumbnail((1000, 760))
            sheet.paste(thumb, (left, top + 55))
    sheet.save(out / "preview/contact_sheet.png")
    return exports


def verify_exports(out, exports):
    require(len(exports) == 8, "Expected four figures in two versions")
    for record in exports:
        require(
            record["minimum_font_pt"] >= (22 if record["mode"] == "poster" else 8),
            "Typography too small",
        )
        require(
            not record["text_overlaps"] and not record["text_outside_canvas"], "Failed layout audit"
        )
        require(
            all((out / "tables" / f"{t}.tsv").is_file() for t in record["source_tables"]),
            "Missing source table",
        )
        probes = json.loads((out / record["value_checks"]).read_text())
        require(
            len(probes) == record["artist_checks"] and all(p["exact_match"] for p in probes),
            "Failed artist validation",
        )
        for filename, checksum in record["files_sha256"].items():
            path = out / filename
            require(sha256(path) == checksum, f"Changed export: {filename}")
            if path.suffix == ".png":
                with Image.open(path) as picture:
                    require(
                        picture.size == tuple(round(v * 300) for v in record["size_inches"]),
                        "PNG dimensions differ",
                    )
                    require(
                        all(abs(v - 300) < 0.1 for v in picture.info.get("dpi", (0, 0))),
                        "PNG is not 300 dpi",
                    )
                    picture.verify()
            elif path.suffix == ".svg":
                ET.parse(path)
            else:
                info = subprocess.check_output(["pdfinfo", str(path)], text=True)
                require(
                    next(v for v in info.splitlines() if v.startswith("Pages:")).split()[-1] == "1",
                    "PDF page count",
                )
