"""Vector data plots from frozen application summaries; no workflow artwork."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib.ticker import PercentFormatter
from PIL import Image, ImageOps

from ._rp1_14_figure_drawing import COLORS, GRAY, INK, Canvas, grid, layout_audit, style

RFU = "#72549A"
INTERVAL_ORDER = {
    "pre_to_last_radiation": 0,
    "last_radiation_to_6weeks": 1,
    "pre_to_6weeks": 2,
}
INTERVAL_LABEL = {
    "pre_to_last_radiation": "P–R",
    "last_radiation_to_6weeks": "R–6w",
    "pre_to_6weeks": "P–6w",
}
TISSUES = {
    "blood": "Blood",
    "bone marrow": "Bone marrow",
    "epithelial lining fluid": "Epithelial lining fluid",
    "inguinal lymph node": "Inguinal lymph node",
    "jejunal epithelium": "Jejunal epithelium",
    "jejunum lamina propria": "Jejunum lamina propria",
    "lung": "Lung",
    "mesenteric lymph node": "Mesenteric lymph node",
    "spleen": "Spleen",
    "thoracic lymph node": "Thoracic lymph node",
}


def representation_key(c, rect):
    ax = c.axes(rect)
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.set_axis_off()
    for x, label, marker, color, face in [
        (0.03, "Receptor TV", "o", GRAY, "white"),
        (0.58, "RFU TV", "D", RFU, RFU),
    ]:
        ax.plot(x, 0.5, marker=marker, color=color, markerfacecolor=face, linestyle="none")
        ax.text(x + 0.045, 0.5, label, va="center", fontsize=c.small)


def paired_vertical(ax, rows, color=INK, open_marker=False):
    for r in rows.itertuples():
        ax.plot([0, 1], [r.d_clone, r.d_group], color=color, alpha=0.35)
        ax.plot(
            [0, 1],
            [r.d_clone, r.d_group],
            linestyle="none",
            marker="s" if open_marker else "o",
            color=color,
            markerfacecolor="white" if open_marker else color,
        )
    for x, col in enumerate(["d_clone", "d_group"]):
        ax.plot([x - 0.17, x + 0.17], [rows[col].median()] * 2, color=INK, linewidth=2.2, zorder=5)
    ax.set(
        xlim=(-0.4, 1.4),
        ylim=(0, 1.055),
        xticks=[0, 1],
        xticklabels=["Receptor", "RFU"],
        yticks=[0, 0.5, 1],
    )
    grid(ax)


def draw_a(c, tables, title=None):
    c.header(
        "A", "GSE190905: paired repertoire change", "Six paired donors · pre/post radiotherapy"
    )
    bottom, height = (0.28, 0.40) if c.compact else (0.57, 0.22)
    for i, subset in enumerate(["all_T", "CD4", "CD8"]):
        ax = c.axes((0.12 + i * 0.30, bottom, 0.20, height))
        part = tables["A_pairs"].query("subset == @subset")
        paired_vertical(ax, part, COLORS.get(subset, INK), subset == "CD8")
        ax.set_title({"all_T": "Total T", "CD4": "CD4", "CD8": "CD8"}[subset] + " · n=6")
        if i == 0:
            ax.set_ylabel("Total variation (TV)")
        else:
            ax.set_yticklabels([])
    if c.compact:
        c.footer(
            "Receptor: nucleotide CDR3/V/J · bars: donor medians\n21,657 / 27,655 cells qualify; ≥100 qualified cells/visit per subset.",
            y=0.025,
        )
        return
    c.text(0.12, 0.465, "Receptor TV (nucleotide CDR3/V/J) · bars: donor medians", color=GRAY)
    coverage = c.axes((0.12, 0.20, 0.29, 0.20))
    for _, part in tables["A_visits"].groupby("donor"):
        coverage.plot(part.time, part.threshold_coverage_of_tcr, "o-", color=GRAY, alpha=0.7)
    coverage.set(
        xlim=(-0.3, 1.3),
        ylim=(0.5, 1),
        xticks=[0, 1],
        xticklabels=["Pre", "Post"],
        yticks=[0.5, 0.75, 1],
        xlabel="Qualified / primary-TRB cells",
        title="Assignment coverage",
    )
    coverage.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    grid(coverage)
    sensitivity = c.axes((0.59, 0.20, 0.35, 0.20))
    conditions = ["cell", "unique", "depth", "dominant_removed"]
    for _, part in tables["A_sensitivity"].groupby("donor"):
        values = part.set_index("condition").reindex(conditions).rfu_tv
        sensitivity.plot(range(4), values, "o-", color=GRAY, alpha=0.55)
    for i, condition in enumerate(conditions):
        source = tables["A_sensitivity"]
        value = source.loc[source.condition.eq(condition), "rfu_tv"].median()
        sensitivity.plot([i - 0.16, i + 0.16], [value] * 2, color=INK, linewidth=2.2)
    sensitivity.set(
        xlim=(-0.4, 3.4),
        ylim=(0.3, 1.03),
        yticks=[0.4, 0.7, 1],
        xticks=range(4),
        xticklabels=["Cell", "Unique\nreceptor", "500\ncells*", "Dominant\nremoved"],
        title="RFU TV sensitivity",
    )
    grid(sensitivity)
    c.footer(
        "21,657 / 27,655 cells qualify · ≥100 qualified cells/visit in all three subsets.\n*50 saved draws/donor; D06 uses 404 cells. Lines pair donors; no new resampling.",
        y=0.025,
    )


def b_order(frame):
    return (
        frame.assign(_interval=frame.interval.map(INTERVAL_ORDER))
        .sort_values(["compartment", "donor", "_interval"])
        .drop(columns="_interval")
        .reset_index(drop=True)
    )


def visit_grid(c, visits, rect, colorbar_rect=None):
    ax = c.axes(rect)
    groups = (
        visits[["donor", "compartment"]].drop_duplicates().sort_values(["donor", "compartment"])
    )
    cmap = plt.get_cmap("Blues")
    for y, r in enumerate(groups.itertuples()):
        for visit in [1, 2, 3]:
            row = visits[
                visits.donor.eq(r.donor)
                & visits.compartment.eq(r.compartment)
                & visits.visit.eq(visit)
            ].iloc[0]
            color = cmap(0.15 + 0.7 * row.coverage) if row.source_available else "#F2F2F2"
            ax.add_patch(
                patches.Rectangle(
                    (visit - 1.5, y - 0.5),
                    1,
                    1,
                    facecolor=color,
                    edgecolor="white",
                    linewidth=2,
                    hatch=None if row.source_available else "///",
                )
            )
            label = (
                f"{int(row.qualified_cells):,}/{int(row.primary_TRB_cells):,}"
                if row.source_available
                else "Missing"
            )
            ax.text(
                visit - 1,
                y,
                label,
                ha="center",
                va="center",
                fontsize=c.small,
                color="white" if row.source_available and row.coverage > 0.65 else INK,
            )
    ax.set(
        xlim=(-0.5, 2.5),
        ylim=(len(groups) - 0.5, -0.5),
        yticks=range(len(groups)),
        yticklabels=[f"{r.donor} {r.compartment}" for r in groups.itertuples()],
        xticks=[0, 1, 2],
        xticklabels=["Pre (P)", "Radiation\nday (R)", "6 weeks\n(6w)"],
    )
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="x", length=0)
    if colorbar_rect:
        bar = c.axes(colorbar_rect)
        # The colorbar has the same restricted palette as the cells.
        from matplotlib.colors import ListedColormap, Normalize

        norm = Normalize(0, 1)
        palette = ListedColormap(cmap(np.linspace(0.15, 0.85, 256)))
        cb = c.fig.colorbar(
            plt.cm.ScalarMappable(norm=norm, cmap=palette),
            cax=bar,
            ticks=[0, 0.5, 1],
            format=PercentFormatter(1, decimals=0),
        )
        cb.outline.set_visible(False)
        cb.solids.set_rasterized(False)


def horizontal_tv(ax, frame, labels=True):
    for y, r in enumerate(frame.itertuples()):
        ax.plot([r.d_group, r.d_clone], [y, y], color="#A8AEB3", linewidth=1.2)
        ax.plot(r.d_clone, y, "o", markerfacecolor="white", color=GRAY)
        ax.plot(r.d_group, y, "D", color=RFU)
    ax.set(
        xlim=(0, 1.045),
        ylim=(len(frame) - 0.5, -0.5),
        xticks=[0, 0.5, 1],
        yticks=range(len(frame)),
        yticklabels=[
            f"{r.donor} {r.compartment} · {INTERVAL_LABEL[r.interval]}" for r in frame.itertuples()
        ]
        if labels
        else [],
        xlabel="Total variation",
    )
    grid(ax, "x")


def draw_b(c, tables, title=None):
    c.header(
        "B",
        "GSE280982: ordered visits and change",
        "11 supported visits · 8 intervals · donors recur across intervals",
    )
    frame = b_order(tables["B_pairs"].query("status == 'analyzed'"))
    if c.compact:
        c.text(0.19, 0.735, "Qualified / primary-TRB cells", color=GRAY)
        visit_grid(c, tables["B_visits"], (0.185, 0.24, 0.36, 0.44))
        ax = c.axes((0.785, 0.24, 0.19, 0.44))
        horizontal_tv(ax, frame)
        representation_key(c, (0.67, 0.705, 0.32, 0.07))
        c.footer(
            "P: pre; R: last radiation day; 6w: six weeks · categorical visit order\nCoverage color: qualified fraction; missing visits are unavailable, never zero.",
            y=0.015,
        )
        return
    c.text(0.24, 0.825, "Qualified / primary-TRB cells; color = qualified fraction", color=GRAY)
    visit_grid(c, tables["B_visits"], (0.24, 0.53, 0.55, 0.265), (0.855, 0.55, 0.023, 0.22))
    tv = c.axes((0.255, 0.14, 0.26, 0.245))
    horizontal_tv(tv, frame)
    representation_key(c, (0.255, 0.403, 0.34, 0.048))
    p = (
        tables["B_persistence"]
        .query("detection_min_cells == 1")
        .set_index(["donor", "compartment", "interval"])
    )
    persistence = c.axes((0.67, 0.14, 0.28, 0.245))
    for i, r in enumerate(frame.itertuples()):
        row = p.loc[(r.donor, r.compartment, r.interval)]
        persistence.plot(
            row.persistent_over_union, i - 0.13, "s", color=GRAY, markerfacecolor="white"
        )
        persistence.plot(row.persistent_without_shared_receptor_fraction, i + 0.13, "D", color=RFU)
    persistence.set(
        xlim=(0, 1),
        ylim=(len(frame) - 0.5, -0.5),
        yticks=[],
        xticks=[0, 0.5, 1],
        xlabel="RFU fraction",
    )
    grid(persistence, "x")
    key = c.axes((0.63, 0.398, 0.34, 0.060))
    key.set_axis_off()
    for y, marker, color, face, label in [
        (0.9, "s", GRAY, "white", "Persistent / union"),
        (0.25, "D", RFU, RFU, "No shared receptor / persistent"),
    ]:
        key.plot(
            0.025,
            y,
            marker=marker,
            color=color,
            markerfacecolor=face,
            linestyle="none",
            transform=key.transAxes,
        )
        key.text(0.075, y, label, transform=key.transAxes, va="center", fontsize=c.small)
    c.footer(
        "Receptor TV: nucleotide CDR3/V/J · persistence: ≥1 cell/endpoint\nRadiation-day tumor depth: 113–157 qualified cells. Persistent RFU ≠ preserved function.",
        y=0.025,
    )


def cross_labels(row):
    if row.dataset == "RP1-14":
        return f"RP1-14 · {row.compartment} · n=6\nReads · 14.95–24.86 years"
    if row.dataset == "GSE190905":
        return "GSE190905 · blood · n=6\nPrimary-TRB cells · pre/post"
    return f"GSE280982 · {row.compartment} · n={row.participants}\nGEX-matched cells · {INTERVAL_LABEL[row.interval]}"


def draw_c(c, tables, title=None):
    c.header(
        "C",
        "Common measurements, distinct sampling units",
        "Study-specific medians · repeated donors · no pooled effect",
    )
    frame = tables["C_summary"]
    bottom, height = (0.20, 0.55) if c.compact else (0.23, 0.49)
    tv = c.axes((0.47, bottom, 0.27, height))
    cancel = c.axes((0.83, bottom, 0.14, height))
    for y, r in enumerate(frame.itertuples()):
        tv.plot([r.receptor_TV_median, r.RFU_TV_median], [y, y], color="#A8AEB3")
        tv.plot(r.receptor_TV_median, y, "o", color=GRAY, markerfacecolor="white")
        tv.plot(r.RFU_TV_median, y, "D", color=RFU)
        cancel.plot(r.cancellation_median, y, "o", color=INK)
        c.text(0.035, bottom + height * (1 - (y + 0.5) / len(frame)), cross_labels(r), va="center")
    for ax in [tv, cancel]:
        ax.set(ylim=(len(frame) - 0.5, -0.5), yticks=[])
        grid(ax, "x")
        for y in [1.5, 2.5]:
            ax.axhline(y, color="#D6DCE1", linewidth=0.6)
    tv.set(xlim=(0, 1.045), xticks=[0, 0.5, 1], xlabel="Total variation")
    cancel.set(xlim=(0, 0.65), xticks=[0, 0.3, 0.6], xlabel="Cancellation")
    representation_key(c, (0.46, 0.755 if c.compact else 0.76, 0.32, 0.06))
    c.text(0.83, 0.831 if c.compact else 0.807, "Aggregation\ncancellation", color=GRAY)
    c.footer(
        "Receptor: nucleotide CDR3/V/J; cancellation = Receptor TV − RFU TV within donor.\nP: pre; R: last radiation day; 6w: six weeks · GEX: gene expression",
        y=0.025,
    )


def draw_d(c, tables, title=None):
    c.header("D", "Regulatory evidence: same-variant overlaps")
    labels = [
        "RFU-QTL\nvariants",
        "Conditional\neQTL overlap",
        "Conditional\ncaQTL overlap",
        "Same-variant\nboth layers",
    ]
    for x, r, label in zip(
        [0.145, 0.385, 0.625, 0.865], tables["D_counts"].itertuples(), labels, strict=True
    ):
        c.text(x, 0.84, f"{r.variants:,}", fontsize=c.number, fontweight="bold", ha="center")
        c.text(x, 0.727, label, ha="center", color=GRAY)
    ax = c.axes((0.045, 0.18, 0.92, 0.39))
    ax.set(xlim=(0, 1), ylim=(6.5, -0.8))
    ax.set_axis_off()
    cols = [
        (0.01, "Variant (chr:pos:ref:alt)"),
        (0.35, "eQTL target"),
        (0.62, "Target class"),
        (0.81, "caQTL\nrecords"),
        (0.96, "RFUs"),
    ]
    for x, label in cols:
        ax.text(
            x,
            -0.75,
            label,
            va="bottom",
            ha="center" if x >= 0.8 else "left",
            fontsize=c.small,
            fontweight="bold",
        )
    for i, row in enumerate(tables["D_matrix"].itertuples()):
        ax.add_patch(
            patches.Rectangle(
                (0, i - 0.43),
                1,
                0.86,
                facecolor="#F4F6F8" if i % 2 == 0 else "white",
                edgecolor="none",
            )
        )
        values = [
            row.variant,
            row.eqtl_target,
            "TCR" if row.target_class == "direct_TCR_gene" else "Non-TCR",
            str(row.caqtl_records),
            str(row.rfu_count),
        ]
        for (x, _), value in zip(cols, values, strict=True):
            ax.text(
                x,
                i,
                value,
                va="center",
                ha="center" if x >= 0.8 else "left",
                fontsize=c.small,
                fontweight="bold" if x == 0.62 else "normal",
            )
    c.footer(
        "eQTL: expression; caQTL: chromatin accessibility · overlaps are not colocalization.\nRFUs recur across rows; caQTL records are distinct target/rank entries. No causal direction.",
        y=0.035,
    )


def draw_e(c, tables, title=None):
    c.header(
        "E",
        "No RFU-added gain in held-out prediction",
        "Six held-out donors · clone-weighted log loss · positive Δ = worse",
    )
    data = tables["E_all_weightings"].query("weighting == 'clone'")
    for i, design in enumerate(["ordinary", "purged"]):
        ax = c.axes((0.12 + i * 0.255, 0.33, 0.165, 0.42))
        part = data[data.design.eq(design)]
        for r in part.itertuples():
            ax.plot([0, 1], [r.baseline_log_loss, r.extended_log_loss], "o-", color=GRAY, alpha=0.7)
        for x, col in enumerate(["baseline_log_loss", "extended_log_loss"]):
            ax.plot([x - 0.18, x + 0.18], [part[col].mean()] * 2, color=INK, linewidth=2.2)
        ax.set(
            xlim=(-0.45, 1.45),
            ylim=(1.7, 2.1),
            yticks=[1.7, 1.8, 1.9, 2.0, 2.1],
            xticks=[0, 1],
            xticklabels=["Baseline", "+RFU"],
            title=design.capitalize(),
        )
        if i == 0:
            ax.set_ylabel("Log loss")
        else:
            ax.set_yticklabels([])
        grid(ax)
    delta = c.axes((0.73, 0.33, 0.245, 0.42))
    for design, offset, marker, color in [
        ("ordinary", -0.12, "o", INK),
        ("purged", 0.12, "D", RFU),
    ]:
        part = data[data.design.eq(design)].sort_values("donor")
        delta.plot(
            part.delta_log_loss,
            np.arange(6) + offset,
            linestyle="none",
            marker=marker,
            color=color,
            label=design.capitalize(),
        )
    delta.axvline(0, color=GRAY, linewidth=1)
    delta.set(
        xlim=(-0.0005, 0.009),
        ylim=(5.5, -0.5),
        xticks=[0, 0.004, 0.008],
        xticklabels=["0", "+.004", "+.008"],
        yticks=range(6),
        yticklabels=part.donor,
        xlabel="Δ log loss\npositive = worse",
    )
    grid(delta, "x")
    delta.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 1.015),
        ncol=2,
        frameon=False,
        handlelength=0.7,
        handletextpad=0.3,
        columnspacing=0.8,
        borderaxespad=0,
    )
    c.text(0.07, 0.195, "Mean Δ: ordinary +0.004216; purged +0.004317", fontweight="bold")
    c.footer(
        "Purged: remove training amino-acid CDR3/V/J identities shared with the test donor.\nBaseline: TRBV, TRBJ and CDR3 length. Bars: equal-donor means; no refitting.",
        y=0.035,
    )


def draw_f(c, tables, title=None):
    c.header(
        "F",
        "Wells: atlas coverage and RFU context",
        "24 atlas donors · 10 tissues · 17 original cell-type labels",
    )
    totals = tables["F_tissues"]
    for x, col, label in [
        (0.17, "atlas_cells", "Atlas cells"),
        (0.49, "primary_trb_cells", "Primary-TRB cells"),
        (0.81, "threshold_cells", "Qualified cells"),
    ]:
        c.text(
            x,
            0.815,
            f"{int(totals[col].sum()):,}",
            fontsize=c.number,
            fontweight="bold",
            ha="center",
        )
        c.text(x, 0.759, label, ha="center", color=GRAY)
    source = tables["F_donor_tissue"]
    tissues, donors = sorted(source.tissue.unique()), sorted(source.donor.unique())
    fractions = source.pivot(
        index="tissue", columns="donor", values="threshold_fraction_of_primary_trb"
    ).reindex(index=tissues, columns=donors)
    ax = c.axes((0.26, 0.445, 0.64, 0.265))
    from matplotlib.colors import Normalize

    norm = Normalize(0, 1)
    cmap = plt.get_cmap("Blues")
    for y, tissue in enumerate(tissues):
        for x, donor in enumerate(donors):
            row = source[source.tissue.eq(tissue) & source.donor.eq(donor)].iloc[0]
            color = (
                cmap(fractions.loc[tissue, donor])
                if row.coverage_state == "primary_TRB_observed"
                else "#E6E8EA"
            )
            ax.add_patch(
                patches.Rectangle(
                    (x - 0.5, y - 0.5), 1, 1, facecolor=color, edgecolor="white", linewidth=0.4
                )
            )
            if row.coverage_state == "atlas_zero_primary_TRB":
                ax.text(x, y, "×", ha="center", va="center", fontsize=c.small)
    ax.set(
        xlim=(-0.5, 23.5),
        ylim=(9.5, -0.5),
        xticks=range(24),
        xticklabels=[d[1:] for d in donors],
        yticks=range(10),
        yticklabels=[TISSUES[t] for t in tissues],
        xlabel="Donor (W01–W24)",
    )
    ax.tick_params(axis="x", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cb = c.fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        cax=c.axes((0.915, 0.465, 0.018, 0.225)),
        ticks=[0, 0.5, 1],
        format=PercentFormatter(1, decimals=0),
    )
    cb.outline.set_visible(False)
    cb.solids.set_rasterized(False)
    c.text(0.26, 0.387, "× atlas sample, zero primary TRB · gray: no atlas sample", color=GRAY)
    c.text(0.26, 0.735, "Qualified / primary-TRB fraction", color=GRAY)
    bar = c.axes((0.26, 0.115, 0.25, 0.215))
    ordered = totals.set_index("tissue").reindex(tissues)
    bar.barh(range(10), ordered.primary_trb_cells / 1000, color="#D2D8DD", height=0.7)
    bar.barh(range(10), ordered.threshold_cells / 1000, color=INK, height=0.35)
    bar.set(
        ylim=(9.5, -0.5),
        yticks=range(10),
        yticklabels=[TISSUES[t] for t in tissues],
        xlabel="Cells (thousands)",
        xticks=[0, 40, 80],
        xlim=(0, 85),
    )
    grid(bar, "x")
    c.text(0.26, 0.35, "Gray: primary; dark: qualified", color=GRAY)
    breadth = c.axes((0.675, 0.115, 0.275, 0.215))
    for dim, marker, color, label in [
        ("tissues", "o", INK, "Tissues"),
        ("cell_types", "s", RFU, "Cell-type labels"),
    ]:
        part = tables["F_breadth_histogram"].query("dimension == @dim")
        breadth.plot(
            part.categories_observed,
            part.rfus,
            marker=marker,
            color=color,
            markerfacecolor="white" if dim == "cell_types" else color,
            label=label,
        )
    breadth.set(
        xlim=(0.5, 17.5),
        xticks=[1, 5, 10, 17],
        xlabel="Categories detected / RFU",
        ylabel="RFUs",
        ylim=(0, 1150),
        yticks=[0, 500, 1000],
    )
    grid(breadth)
    breadth.legend(
        loc="lower left",
        bbox_to_anchor=(-0.03, 1.04),
        frameon=False,
        ncol=1,
        handlelength=1,
        borderaxespad=0,
        labelspacing=0.2,
    )
    c.footer(
        "4,928 RFUs; 1,471 meet ≥50 cells across ≥4 donors. Detection breadth ≠ specificity.",
        y=0.025,
    )


DRAW = {"A": draw_a, "B": draw_b, "C": draw_c, "D": draw_d, "E": draw_e, "F": draw_f}


def export(fig, out, stem, key, mode, config, title, **extra):
    settings = config["exports"]
    audit = layout_audit(fig, settings[f"{mode}_min_font_pt"])
    files = []
    for ext in ["pdf", "svg", "png"]:
        path = out / f"{stem}.{ext}"
        metadata = {"Creator": "scRFU completed-application figure pipeline", "Title": title}
        if ext == "pdf":
            metadata.update(CreationDate=None, ModDate=None)
        elif ext == "svg":
            metadata["Date"] = None
        fig.savefig(path, format=ext, dpi=settings["dpi"], metadata=metadata)
        files.append(str(path.relative_to(out)))
    preview = out / "preview" / (Path(stem).name + ".png")
    fig.savefig(preview, dpi=settings["preview_dpi"])
    with Image.open(preview) as img:
        ImageOps.grayscale(img.convert("RGB")).save(
            preview.with_name(preview.stem + "_grayscale.png")
        )
    size = fig.get_size_inches().tolist()
    plt.close(fig)
    return {
        "figure": key,
        "mode": mode,
        "files": files,
        "preview": str(preview.relative_to(out)),
        "grayscale_preview": str(
            preview.with_name(preview.stem + "_grayscale.png").relative_to(out)
        ),
        "size_inches": size,
        "size_mm": [round(v * 25.4, 2) for v in size],
        "dpi": settings["dpi"],
        "vector_pdf_svg": True,
        **audit,
        **extra,
    }


def contact_sheet(out, mode):
    tiles = []
    for key in "ABCDEF":
        with Image.open(out / "preview" / f"cross_application_{key}_{mode}.png") as img:
            tile = img.convert("RGB")
            tile.thumbnail((820, 680))
            tiles.append(tile)
    sheet = Image.new("RGB", (1680, 2100), "#DDE3E7")
    for i, tile in enumerate(tiles):
        sheet.paste(tile, (12 + (i % 2) * 840, 12 + (i // 2) * 700))
    sheet.save(out / "preview" / f"contact_sheet_{mode}.png")


def render_all(out, tables, config, captions):
    exports = []
    for mode in ["manuscript", "poster"]:
        with plt.rc_context({**style(mode), "svg.hashsalt": "cross-application-reusable-v1"}):
            for key, draw in DRAW.items():
                fig = plt.figure(figsize=config["exports"][f"{mode}_sizes_inches"][key])
                draw(Canvas(fig, (0, 0, 1, 1), mode), tables)
                exports.append(
                    export(
                        fig,
                        out,
                        f"figures/{mode}/cross_application_{key}_{mode}",
                        key,
                        mode,
                        config,
                        captions[key]["title"],
                    )
                )
            contact_sheet(out, mode)
    with plt.rc_context({**style("manuscript"), "svg.hashsalt": "cross-application-reusable-v1"}):
        for key, size in [("A", (7.5, 2.8)), ("B", (7.5, 3.4))]:
            fig = plt.figure(figsize=size)
            DRAW[key](Canvas(fig, (0, 0, 1, 1), "manuscript", compact=True), tables)
            exports.append(
                export(
                    fig,
                    out,
                    f"figures/panels/cross_application_{key}_compact",
                    key,
                    "manuscript",
                    config,
                    captions[key]["title"],
                    variant="compact",
                )
            )
        fig = plt.figure(figsize=config["exports"]["manuscript_assembly_inches"])
        fig.text(
            0.035,
            0.994,
            captions["manuscript_assembly"]["title"],
            fontsize=12,
            fontweight="bold",
            va="top",
        )
        for key, rect in [
            ("A", (0.01, 0.696, 0.98, 0.266)),
            ("B", (0.01, 0.367, 0.98, 0.323)),
            ("C", (0.01, 0.004, 0.98, 0.354)),
        ]:
            DRAW[key](Canvas(fig, rect, "manuscript", compact=True), tables)
        exports.append(
            export(
                fig,
                out,
                "figures/manuscript/cross_application_manuscript_assembly",
                "manuscript_assembly",
                "manuscript",
                config,
                captions["manuscript_assembly"]["title"],
            )
        )
    with plt.rc_context({**style("poster"), "svg.hashsalt": "cross-application-reusable-v1"}):
        fig = plt.figure(figsize=config["exports"]["poster_block_inches"])
        fig.text(
            0.03, 0.992, captions["poster_block"]["title"], fontsize=34, fontweight="bold", va="top"
        )
        for key, rect in [
            ("A", (0.01, 0.49, 0.48, 0.455)),
            ("B", (0.515, 0.49, 0.48, 0.455)),
            ("D", (0.01, 0.045, 0.59, 0.375)),
        ]:
            DRAW[key](Canvas(fig, rect, "poster"), tables)
        fig.text(
            0.03,
            0.012,
            "RFU: receptor functional unit · TV: total variation · TRB: T-cell receptor beta",
            fontsize=18,
        )
        exports.append(
            export(
                fig,
                out,
                "figures/poster/cross_application_poster_block",
                "poster_block",
                "poster",
                config,
                captions["poster_block"]["title"],
                biorender_blank_rectangle_fraction=config["exports"][
                    "poster_biorender_blank_rectangle_fraction"
                ],
            )
        )
    return exports
