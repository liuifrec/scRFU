"""Matplotlib-only vector panels shared by RP1-14 manuscript/poster exports."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib.text import Text
from matplotlib.ticker import PercentFormatter
from matplotlib.transforms import Bbox
from PIL import Image, ImageOps

COLORS = {"CD4": "#0072B2", "CD8": "#C15A00"}
MARKERS = {"CD4": "o", "CD8": "s"}
INK = "#18232D"
GRAY = "#63717B"
PALE = "#F3F6F8"
CONDITIONS = ["read_weighted", "unique_receptor", "500_reads"]


def style(mode: str) -> dict:
    poster = mode == "poster"
    return {
        "font.family": "DejaVu Sans",
        "font.size": 20 if poster else 9,
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.edgecolor": GRAY,
        "axes.labelsize": 20 if poster else 8,
        "axes.titlesize": 22 if poster else 9,
        "xtick.labelsize": 18 if poster else 8,
        "ytick.labelsize": 18 if poster else 8,
        "legend.fontsize": 18 if poster else 8,
        "axes.linewidth": 1.1 if poster else 0.6,
        "lines.linewidth": 1.5 if poster else 0.8,
        "lines.markersize": 7 if poster else 3.8,
        "xtick.major.size": 4 if poster else 2.5,
        "ytick.major.size": 0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.axisbelow": True,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "path",
        "svg.hashsalt": "rp1-14-reusable-v1",
        "figure.max_open_warning": 0,
    }


@dataclass
class Canvas:
    fig: plt.Figure
    rect: tuple[float, float, float, float]
    mode: str
    compact: bool = False

    @property
    def poster(self) -> bool:
        return self.mode == "poster"

    @property
    def small(self) -> int:
        return 18 if self.poster else 8

    @property
    def normal(self) -> int:
        return 20 if self.poster else 9

    @property
    def number(self) -> int:
        return 28 if self.poster else (11 if self.compact else 14)

    def point(self, x: float, y: float) -> tuple[float, float]:
        left, bottom, width, height = self.rect
        return left + x * width, bottom + y * height

    def text(self, x: float, y: float, s: str, **kwargs) -> Text:
        defaults = {"fontsize": self.small, "va": "top", "ha": "left"}
        defaults.update(kwargs)
        return self.fig.text(*self.point(x, y), s, **defaults)

    def axes(self, rect: tuple[float, float, float, float]) -> plt.Axes:
        left, bottom = self.point(rect[0], rect[1])
        return self.fig.add_axes([left, bottom, rect[2] * self.rect[2], rect[3] * self.rect[3]])

    def header(self, letter: str, title: str, subtitle: str = "") -> None:
        self.text(0.022, 0.973, letter, fontsize=26 if self.poster else 12, fontweight="bold")
        self.text(0.07, 0.968, title, fontsize=26 if self.poster else 11, fontweight="bold")
        if subtitle:
            self.text(0.07, 0.873, subtitle, color=GRAY)

    def footer(self, text: str, y: float = 0.026, **kwargs) -> None:
        self.text(0.07, y, text, va="bottom", **kwargs)


def compartment_key(c: Canvas, rect=(0.1, 0.78, 0.65, 0.06), *, controls=False) -> None:
    ax = c.axes(rect)
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.set_axis_off()
    if controls:
        ax.text(0, 0.5, "RFU:", fontsize=c.small, va="center")
    for comp, x in [("CD4", 0.11 if controls else 0.025), ("CD8", 0.28)]:
        ax.plot(
            x,
            0.5,
            marker=MARKERS[comp],
            color=COLORS[comp],
            markerfacecolor=COLORS[comp] if comp == "CD4" else "white",
            linestyle="none",
        )
        ax.text(x + 0.035, 0.5, comp, fontsize=c.small, va="center", color=COLORS[comp])
    if controls:
        ax.plot(0.46, 0.5, marker="D", markerfacecolor="white", color=GRAY, linestyle="none")
        ax.text(0.495, 0.5, "Matched maps: median / range", va="center", fontsize=c.small)


def points(ax, x, y, comp: str, *, size=None, alpha=1, zorder=4) -> None:
    ax.scatter(
        x,
        y,
        marker=MARKERS[comp],
        s=size or (plt.rcParams["lines.markersize"] ** 2),
        facecolor=COLORS[comp] if comp == "CD4" else "white",
        edgecolor=COLORS[comp],
        linewidth=1,
        alpha=alpha,
        zorder=zorder,
    )


def grid(ax: plt.Axes, axis="y") -> None:
    ax.grid(axis=axis, color="#E4E9EC", linewidth=0.6)


def draw_a(c: Canvas, tables: dict, title: str) -> None:
    totals = tables["A_totals"].set_index("metric").value
    span = f"{totals.follow_up_years_min:.2f}–{totals.follow_up_years_max:.2f} years"
    c.header("A", title, f"CD4 and CD8 at every visit · 36 samples · follow-up {span}")
    samples, depth = tables["A_samples"], tables["A_donor_depth"]
    donors = sorted(samples.donor.unique())
    ax = c.axes((0.10, 0.35, 0.48, 0.39 if c.compact else 0.42))
    for i, donor in enumerate(donors):
        for comp, offset in [("CD4", -0.12), ("CD8", 0.12)]:
            part = samples[samples.donor.eq(donor) & samples.compartment.eq(comp)]
            ax.plot(
                part.elapsed_years, np.full(len(part), i + offset), color=COLORS[comp], alpha=0.45
            )
            points(ax, part.elapsed_years, np.full(len(part), i + offset), comp)
    ax.set(
        yticks=range(6),
        yticklabels=donors,
        ylim=(5.6, -0.6),
        xlim=(-0.8, 26),
        xticks=[0, 5, 10, 15, 20, 25],
        xlabel="Years from first visit",
    )
    ax.spines["left"].set_visible(False)
    compartment_key(c, (0.11, 0.774, 0.40, 0.055))
    bar = c.axes((0.70, 0.35, 0.265, 0.39 if c.compact else 0.42))
    for i, donor in enumerate(donors):
        for comp, offset in [("CD4", -0.17), ("CD8", 0.17)]:
            row = depth[depth.donor.eq(donor) & depth.compartment.eq(comp)].iloc[0]
            bar.barh(i + offset, row.all_source_reads / 1e6, height=0.27, color="#DDE3E7", zorder=2)
            bar.barh(
                i + offset,
                row.fixed_primary_reads / 1e6,
                height=0.27,
                facecolor=COLORS[comp] if comp == "CD4" else "white",
                edgecolor=COLORS[comp],
                hatch=None if comp == "CD4" else "///",
                linewidth=0.65,
                zorder=3,
            )
    bar.set(ylim=(5.6, -0.6), yticks=[], xlabel="Reads / donor (millions)")
    bar.set_xlim(0, np.ceil(depth.all_source_reads.max() / 1e6))
    bar.xaxis.set_major_locator(plt.MaxNLocator(3, integer=True))
    bar.spines["left"].set_visible(False)
    c.text(0.70, 0.81, "Gray: source; color: qualified", fontsize=c.small)
    metrics = [
        ("source_reads", "Source reads"),
        ("primary_mapped_reads", "Qualified reads"),
        ("primary_unique_receptors", "Distinct nucleotide\nCDR3/V/J identities"),
        ("primary_unique_rfus", "Receptor functional\nunits (RFUs)"),
    ]
    for x, (key, label) in zip([0.14, 0.39, 0.65, 0.90], metrics, strict=True):
        c.text(x, 0.19, f"{int(totals[key]):,}", ha="center", fontsize=c.number, fontweight="bold")
        c.text(x, 0.105 if c.compact else 0.12, label, ha="center", color=GRAY)
    if not c.compact:
        c.footer(
            "Qualified: fixed historical dictionary, score ≥0.6. Sequencing reads are not cell counts."
        )


def representation_axis(ax, pairs, comp: str, *, poster: bool) -> None:
    metrics = ["receptor_tv", "rfu_tv"] if poster else ["receptor_tv", "rfu_tv", "trbv_tv", "vj_tv"]
    labels = (
        ["Receptor TV\n(nucleotide\nCDR3/V/J)", "RFU TV"]
        if poster
        else ["Receptor TV\n(nucleotide\nCDR3/V/J)", "RFU TV", "TRBV TV", "V/J TV"]
    )
    x = np.arange(len(metrics))
    group = pairs[pairs.compartment.eq(comp)].sort_values("donor")
    for i, row in enumerate(group.to_dict(orient="records")):
        values = [row[metric] for metric in metrics]
        xx = x + (i - 2.5) * 0.025
        ax.plot(xx, values, color=COLORS[comp], alpha=0.36)
        points(ax, xx, values, comp, alpha=0.9)
    for i, metric in enumerate(metrics):
        ax.plot(
            [i - 0.17, i + 0.17],
            [group[metric].median()] * 2,
            color=INK,
            lw=2.8 if poster else 1.8,
            zorder=6,
        )
    ax.set(
        ylim=(0, 1.02),
        xlim=(-0.35, len(metrics) - 0.65),
        yticks=[0, 0.5, 1],
        xticks=x,
        xticklabels=labels,
    )
    ax.set_title(f"{comp} · 6 donors", loc="left", color=COLORS[comp], fontweight="bold", pad=5)
    grid(ax)


def draw_b(c: Canvas, tables: dict, title: str) -> None:
    c.header(
        "B", title, "Earliest → latest · read weighted · lines: same donor; black bars: medians"
    )
    for i, comp in enumerate(COLORS):
        ax = c.axes(
            (0.105 + i * 0.47, 0.32 if c.compact else 0.27, 0.375, 0.415 if c.compact else 0.465)
        )
        representation_axis(ax, tables["B_pairs"], comp, poster=c.poster)
        if i == 0:
            ax.set_ylabel("Total variation (TV)")
    if c.poster:
        summary = tables["B_pairs"].groupby("compartment")[["receptor_tv", "rfu_tv"]].median()
        c.footer(
            "Donor medians: "
            + "    ".join(
                f"{comp}  {row.receptor_tv:.3f} → {row.rfu_tv:.3f}"
                for comp, row in summary.iterrows()
            ),
            y=0.11,
            fontweight="bold",
        )
    c.footer(
        "Lower group TV partly reflects aggregation."
        if c.compact
        else "TV: 0 = identical, 1 = disjoint. Lower group TV partly reflects aggregation."
    )


def controls_axis(ax, frame, comp: str) -> None:
    part = frame[frame.compartment.eq(comp)].sort_values("donor")
    for i, row in enumerate(part.itertuples()):
        ax.plot(
            [row.observed_cancellation, row.control_median],
            [i - 0.11, i + 0.11],
            color="#B1BCC4",
            lw=0.8,
        )
        ax.plot([row.control_q025, row.control_q975], [i + 0.11] * 2, color=COLORS[comp], lw=1.4)
        ax.plot(
            row.control_median,
            i + 0.11,
            marker="D",
            markerfacecolor="white",
            color=COLORS[comp],
            markersize=plt.rcParams["lines.markersize"] * 0.8,
        )
        points(ax, [row.observed_cancellation], [i - 0.11], comp)
    ax.set(
        yticks=range(len(part)),
        yticklabels=part.donor,
        ylim=(5.6, -0.6),
        xlim=(0, 0.7),
        xticks=[0, 0.2, 0.4, 0.6],
        xlabel="Aggregation cancellation",
    )
    ax.set_title(f"{comp} · 6 donors", loc="left", color=COLORS[comp], fontweight="bold", pad=5)
    ax.spines["left"].set_visible(False)
    grid(ax, "x")


def draw_c(c: Canvas, tables: dict, title: str) -> None:
    c.header("C", title, "Aggregation cancellation = Receptor TV (nucleotide CDR3/V/J) − RFU TV")
    compartment_key(c, (0.10, 0.748, 0.84, 0.055), controls=True)
    for i, comp in enumerate(COLORS):
        ax = c.axes(
            (0.105 + i * 0.48, 0.335 if c.compact else 0.255, 0.37, 0.35 if c.compact else 0.425)
        )
        controls_axis(ax, tables["C_donors"], comp)
    c.footer("Similar cancellation is not uniquely RFU-specific.", y=0.085, fontweight="bold")
    c.footer(
        "Ranges: 2.5–97.5% of 30 matched maps; technical variation, not donor uncertainty.", y=0.015
    )


def persistence_axis(ax, frame, metric: str, c: Canvas) -> None:
    pivot = frame.pivot(index="donor", columns="compartment", values=metric)
    for i, (_, row) in enumerate(pivot.iterrows()):
        offset = (i - 2.5) * 0.032
        ax.plot(np.array([0, 1]) + offset, [row.CD4, row.CD8], color="#AAB5BD", alpha=0.7)
        for x, comp in enumerate(COLORS):
            points(ax, [x + offset], [row[comp]], comp)
    for x, comp in enumerate(COLORS):
        median = pivot[comp].median()
        ax.plot(
            [x - 0.19, x + 0.19], [median] * 2, color=INK, lw=2.8 if c.poster else 1.8, zorder=6
        )
        ax.text(
            x,
            0.96,
            f"{100 * median:.1f}%",
            ha="center",
            va="top",
            fontsize=c.number,
            color=COLORS[comp],
            fontweight="bold",
        )
    ax.set(
        xlim=(-0.35, 1.35),
        ylim=(0, 1.02),
        xticks=[0, 1],
        xticklabels=["CD4", "CD8"],
        yticks=[0, 0.25, 0.5, 0.75, 1],
    )
    ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    grid(ax)


def schematic_axis(ax: plt.Axes, frame, c: Canvas) -> None:
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.set_axis_off()
    ax.text(0.5, 0.98, "Explanatory schematic", ha="center", va="top", fontsize=c.small, color=GRAY)
    ax.text(
        0.5,
        0.81,
        "One RFU across two visits",
        ha="center",
        va="center",
        fontsize=c.normal,
        fontweight="bold",
    )
    for visit, x, label in [("early", 0.05, "Early"), ("late", 0.57, "Late")]:
        ax.add_patch(
            patches.FancyBboxPatch(
                (x, 0.35),
                0.37,
                0.32,
                boxstyle="round,pad=0.015,rounding_size=0.03",
                facecolor=PALE,
                edgecolor=GRAY,
                lw=1.2,
            )
        )
        ax.text(x + 0.185, 0.71, label, ha="center", va="center", fontsize=c.small)
        letters = frame[frame.visit.eq(visit)].arbitrary_receptor_label.tolist()
        for pos, letter in zip([x + 0.09, x + 0.28], letters, strict=True):
            ax.text(
                pos,
                0.51,
                letter,
                ha="center",
                va="center",
                fontsize=c.normal,
                fontweight="bold",
                bbox={"boxstyle": "circle,pad=0.28", "facecolor": "white", "edgecolor": GRAY},
            )
    ax.annotate(
        "",
        xy=(0.545, 0.51),
        xytext=(0.445, 0.51),
        arrowprops={"arrowstyle": "->", "color": GRAY, "lw": 1.5},
    )
    ax.text(0.5, 0.20, "Different observed receptors", ha="center", va="center", fontsize=c.small)
    ax.text(
        0.5,
        0.06,
        "Arbitrary labels; not measured data",
        ha="center",
        va="center",
        fontsize=c.small,
        color=GRAY,
    )


def draw_d(c: Canvas, tables: dict, title: str) -> None:
    c.header("D", title, "Detection ≥5 qualified reads in both endpoint visits · 6 paired donors")
    if c.poster or c.compact:
        ax = c.axes((0.115, 0.215, 0.33, 0.545))
        persistence_axis(ax, tables["D_donors"], "fraction_persistent_without_shared_clone", c)
        ax.set_ylabel("Persistent RFUs without a\nshared observed receptor")
        schematic_axis(c.axes((0.515, 0.195, 0.46, 0.57)), tables["D_schematic"], c)
    else:
        for i, (metric, label) in enumerate(
            [
                ("persistent_fraction_union", "Persistent / union detected RFUs"),
                (
                    "fraction_persistent_without_shared_clone",
                    "No shared receptor / persistent RFUs",
                ),
            ]
        ):
            ax = c.axes((0.105 + i * 0.49, 0.235, 0.36, 0.51))
            persistence_axis(ax, tables["D_donors"], metric, c)
            ax.set_title(label, loc="center", fontsize=c.small, pad=9)
            if i == 0:
                ax.set_ylabel("Fraction of RFUs")
    c.footer(
        "Observed persistence does not establish preserved antigen function.",
        y=0.036,
        fontweight="bold",
    )
    if not c.compact:
        c.footer("Lines pair donors; black bars and percentages show donor medians.", y=0.102)


def sensitivity_axis(ax, frame, metric: str, comp: str, c: Canvas, *, intervals: bool) -> None:
    part = frame[frame.compartment.eq(comp)]
    for i, (_, group) in enumerate(part.groupby("donor")):
        rows = group.set_index("condition").loc[CONDITIONS]
        xx = np.arange(3) + (i - 2.5) * 0.027
        ax.plot(xx, rows[metric], color=COLORS[comp], alpha=0.32)
        points(ax, xx, rows[metric], comp, alpha=0.7, size=26 if c.poster else 11)
        if intervals:
            ax.plot(
                [xx[-1]] * 2,
                [rows.loc["500_reads", f"{metric}_lo"], rows.loc["500_reads", f"{metric}_hi"]],
                color=COLORS[comp],
                alpha=0.55,
            )
    meds = part.groupby("condition")[metric].median().reindex(CONDITIONS)
    if c.poster:
        ax.plot(
            range(3), meds, color=COLORS[comp], lw=3, ls="-" if comp == "CD4" else "--", zorder=5
        )
        points(ax, range(3), meds, comp, size=100, zorder=6)
    else:
        for x, med in enumerate(meds):
            ax.plot([x - 0.14, x + 0.14], [med] * 2, color=INK, lw=1.6, zorder=6)
    ax.set(
        ylim=(0, 1.02),
        yticks=[0, 0.5, 1],
        xlim=(-0.3, 2.3),
        xticks=range(3),
        xticklabels=["Read\nweighted", "Unique\nreceptor", "500 reads\nper visit"],
    )
    grid(ax)


def draw_e(c: Canvas, tables: dict, title: str) -> None:
    if c.poster:
        narrow = c.rect[2] * c.fig.get_figwidth() < 11
        c.header(
            "E",
            "Sampling and weighting\nchange the estimate" if narrow else title,
            "RFU TV · the same six donors" if not narrow else "",
        )
        if narrow:
            c.text(0.07, 0.83, "RFU TV · the same six donors")
        compartment_key(c, (0.18, 0.755, 0.65, 0.06))
        ax = c.axes((0.14, 0.275, 0.82, 0.44))
        for comp in COLORS:
            sensitivity_axis(ax, tables["E_donors"], "rfu_tv", comp, c, intervals=False)
        ax.set_ylabel("RFU TV")
        c.footer("Thin lines: donors; heavy lines: medians.", y=0.14)
        c.footer(
            "500-read points: within-donor medians\nof 50 frozen draws; not donor uncertainty.",
            y=0.035,
        )
    else:
        c.header(
            "E",
            title,
            "Earliest → latest · same qualified receptor universe · 6 donors per compartment",
        )
        for row, (metric, label) in enumerate(
            [("rfu_tv", "RFU TV"), ("cancellation", "Aggregation\ncancellation")]
        ):
            for i, comp in enumerate(COLORS):
                ax = c.axes((0.10 + i * 0.49, 0.52 - row * 0.325, 0.36, 0.24))
                sensitivity_axis(ax, tables["E_donors"], metric, comp, c, intervals=True)
                if i == 0:
                    ax.set_ylabel(label)
                if row == 0:
                    ax.set_title(comp, color=COLORS[comp], loc="left", fontweight="bold", pad=4)
                    ax.set_xticklabels([])
        c.footer(
            "500-read whiskers: 2.5–97.5% of 50 frozen draws; technical variation, not donor uncertainty.",
            y=0.025,
        )


def draw_f(c: Canvas, tables: dict, title: str) -> None:
    checks = tables["F_checks"].set_index("check").value
    c.header("F", title, "Completed historical row-alignment audit · no assignment or assay rerun")
    cards = [
        (
            f"{int(checks.vector_entries_unchanged):,}",
            "RFU / score entries unchanged",
            f"{int(checks.amino_acid_labels_repaired):,} amino-acid labels repaired",
        ),
        (
            f"{int(checks.rfu_conflicts_before):,} → {int(checks.rfu_conflicts_after)}",
            "Repeated-amino-acid RFU conflicts",
            f"Score conflicts: {int(checks.score_conflicts_before):,} → {int(checks.score_conflicts_after)}",
        ),
        (
            f"{int(checks.rfu_label_matches):,} / 1,600",
            "Label and threshold parity matches",
            f"{int(checks.near_threshold_assays):,} scores in [0.59, 0.61) included",
        ),
        (
            rf"${checks.max_abs_score_difference / 1e-15:.3f}\times10^{{-15}}$",
            "Maximum absolute score discrepancy",
            f"Frozen tolerance: {checks.score_tolerance:.0e}",
        ),
    ]
    for i, (value, label, detail) in enumerate(cards):
        left, bottom = 0.055 + (i % 2) * 0.48, 0.475 - (i // 2) * 0.305
        ax = c.axes((left, bottom, 0.43, 0.255))
        ax.set_axis_off()
        ax.add_patch(
            patches.FancyBboxPatch(
                (0, 0),
                1,
                1,
                transform=ax.transAxes,
                boxstyle="round,pad=0.0,rounding_size=0.04",
                facecolor=PALE,
                edgecolor="#D9E1E6",
            )
        )
        ax.text(0.055, 0.78, value, fontsize=34 if c.poster else 19, fontweight="bold", va="center")
        ax.text(0.055, 0.43, label, fontsize=c.small, fontweight="bold", va="center")
        ax.text(0.055, 0.18, detail, fontsize=c.small, color=GRAY, va="center")
    c.footer(
        "Bounded repair evidence; the historical invocation manifest remains unavailable.", y=0.04
    )


DRAW = {"A": draw_a, "B": draw_b, "C": draw_c, "D": draw_d, "E": draw_e, "F": draw_f}


def layout_audit(fig, minimum: float) -> dict:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bounds = fig.bbox
    outside, sizes, boxes = [], [], []
    artists = list(fig.texts)
    for ax in fig.axes:
        artists.extend(ax.texts)
        artists.extend([ax.title, ax._left_title, ax._right_title])
        if ax.axison:
            for axis in [ax.xaxis, ax.yaxis]:
                if not axis.get_visible():
                    continue
                artists.extend([axis.label, axis.get_offset_text()])
                low, high = sorted(axis.get_view_interval())
                # Locators also create unused ticks beyond the view interval.
                for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                    if low - 1e-9 <= tick.get_loc() <= high + 1e-9:
                        artists.extend([tick.label1, tick.label2])
    for artist in artists:
        if not artist.get_visible() or not artist.get_text().strip():
            continue
        box = artist.get_window_extent(renderer)
        boxes.append((artist.get_text(), box))
        sizes.append(artist.get_fontsize())
        if (
            box.x0 < bounds.x0 - 1
            or box.y0 < bounds.y0 - 1
            or box.x1 > bounds.x1 + 1
            or box.y1 > bounds.y1 + 1
        ):
            outside.append(artist.get_text())
    if sizes and min(sizes) < minimum:
        raise ValueError(f"Figure font below {minimum} pt: {min(sizes)}")
    if outside:
        raise ValueError(f"Text outside figure canvas: {outside}")
    overlaps = []
    for (a, first), (b, second) in combinations(boxes, 2):
        intersection = Bbox.intersection(first, second)
        if (
            intersection is not None
            and intersection.width > 1
            and intersection.height > 1
            and intersection.width
            * intersection.height
            / min(first.width * first.height, second.width * second.height)
            > 0.08
        ):
            overlaps.append([a, b])
    if overlaps:
        raise ValueError(f"Overlapping text labels: {overlaps}")
    return {
        "minimum_font_pt": min(sizes),
        "required_minimum_font_pt": minimum,
        "text_outside_canvas": outside,
        "text_overlaps": overlaps,
    }


def export(fig, out: Path, stem: str, key: str, mode: str, config: dict, *, title: str) -> dict:
    settings = config["exports"]
    audit = layout_audit(fig, settings[f"{mode}_min_font_pt"])
    paths = []
    for extension in settings["formats"]:
        path = out / f"{stem}.{extension}"
        metadata = {"Creator": "scRFU reusable RP1-14 figure pipeline", "Title": title}
        if extension == "pdf":
            metadata.update({"CreationDate": None, "ModDate": None})
        elif extension == "svg":
            metadata["Date"] = None
        fig.savefig(path, format=extension, dpi=settings["png_dpi"], metadata=metadata)
        paths.append(str(path.relative_to(out)))
    preview = out / "preview" / (Path(stem).name + ".png")
    fig.savefig(preview, dpi=settings["preview_dpi"])
    with Image.open(preview) as img:
        gray = ImageOps.grayscale(img.convert("RGB"))
        gray.save(preview.with_name(preview.stem + "_grayscale.png"))
    size = fig.get_size_inches().tolist()
    plt.close(fig)
    return {
        "figure": key,
        "mode": mode,
        "files": paths,
        "preview": str(preview.relative_to(out)),
        "grayscale_preview": str(
            preview.with_name(preview.stem + "_grayscale.png").relative_to(out)
        ),
        "size_inches": size,
        "size_mm": [round(x * 25.4, 2) for x in size],
        "dpi": settings["png_dpi"],
        "vector_pdf_svg": True,
        **audit,
    }


def contact_sheet(out: Path, mode: str) -> None:
    images = []
    for key in "ABCDEF":
        with Image.open(out / "preview" / f"rp1_14_{key}_{mode}.png") as img:
            copy = img.convert("RGB")
            copy.thumbnail((850, 565))
            images.append(copy)
    width = max(img.width for img in images)
    height = max(img.height for img in images)
    sheet = Image.new("RGB", (2 * width + 36, 3 * height + 48), "#DAE0E4")
    for i, img in enumerate(images):
        sheet.paste(img, (12 + (i % 2) * (width + 12), 12 + (i // 2) * (height + 12)))
    sheet.save(out / "preview" / f"contact_sheet_{mode}.png")


def render_all(out: Path, tables: dict, config: dict, captions: dict) -> list[dict]:
    exports = []
    for mode in ["manuscript", "poster"]:
        with plt.rc_context(style(mode)):
            for key, draw in DRAW.items():
                fig = plt.figure(figsize=config["exports"][f"{mode}_inches"])
                draw(Canvas(fig, (0, 0, 1, 1), mode), tables, captions[key]["title"])
                stem = f"figures/{mode}/rp1_14_{key}_{mode}"
                exports.append(
                    export(fig, out, stem, key, mode, config, title=captions[key]["title"])
                )
            fig = plt.figure(figsize=(7, 3.2) if mode == "manuscript" else (12, 5.2))
            c = Canvas(fig, (0, 0, 1, 1), mode)
            schematic_axis(c.axes((0.04, 0.02, 0.92, 0.94)), tables["D_schematic"], c)
            exports.append(
                export(
                    fig,
                    out,
                    f"figures/panels/rp1_14_turnover_schematic_{mode}",
                    "D",
                    mode,
                    config,
                    title="Explanatory schematic; not measured data",
                )
            )
            fig = plt.figure(figsize=(4, 0.7) if mode == "manuscript" else (8, 1.3))
            c = Canvas(fig, (0, 0, 1, 1), mode)
            compartment_key(c, (0.05, 0.1, 0.90, 0.8))
            exports.append(
                export(
                    fig,
                    out,
                    f"figures/panels/rp1_14_compartment_key_{mode}",
                    "A",
                    mode,
                    config,
                    title="Reusable CD4/CD8 encoding",
                )
            )
            contact_sheet(out, mode)
    with plt.rc_context(style("manuscript")):
        fig = plt.figure(figsize=config["exports"]["manuscript_assembly_inches"])
        fig.text(
            0.035,
            0.995,
            captions["manuscript_main"]["title"],
            fontsize=12,
            fontweight="bold",
            va="top",
        )
        for key, bottom in zip("ABCD", [0.737, 0.497, 0.253, 0.009], strict=True):
            DRAW[key](
                Canvas(fig, (0.01, bottom, 0.98, 0.225), "manuscript", compact=True),
                tables,
                captions[key]["title"],
            )
        exports.append(
            export(
                fig,
                out,
                "figures/manuscript/rp1_14_manuscript_main",
                "manuscript_main",
                "manuscript",
                config,
                title=captions["manuscript_main"]["title"],
            )
        )
    with plt.rc_context(style("poster")):
        fig = plt.figure(figsize=config["exports"]["poster_assembly_inches"])
        fig.text(
            0.03, 0.994, captions["poster_block"]["title"], fontsize=34, fontweight="bold", va="top"
        )
        positions = {
            "A": (0.01, 0.525, 0.48, 0.425),
            "C": (0.515, 0.525, 0.48, 0.425),
            "D": (0.01, 0.055, 0.62, 0.425),
            "E": (0.665, 0.055, 0.32, 0.425),
        }
        for key, rect in positions.items():
            DRAW[key](Canvas(fig, rect, "poster"), tables, captions[key]["title"])
        fig.text(
            0.03,
            0.012,
            "RFU: receptor functional unit   ·   TV: total variation   ·   The same six donors contribute CD4 and CD8",
            fontsize=18,
        )
        exports.append(
            export(
                fig,
                out,
                "figures/poster/rp1_14_poster_block",
                "poster_block",
                "poster",
                config,
                title=captions["poster_block"]["title"],
            )
        )
    return exports
