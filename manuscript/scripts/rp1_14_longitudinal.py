"""Published-only RP1-14 longitudinal measurements using repaired historical RFUs.

Apply a fixed, recovered amino-acid-to-RFU dictionary to *all* observed productive
receptor counts at every visit. This avoids treating a receptor below one visit's
historical top-10000 cutoff as absent. Unmapped read mass remains explicit. No
receptors are newly assigned and no previous public analyses are rerun.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import inspect
import itertools
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import scrfu
from scrfu.io import file_sha256
from scrfu.tl import longitudinal_similarity, multiscale_repertoire_change, permute_fixed_groups

from ._figure_common import clean_axis, panel_label, save_figure, set_style
from .radiation_methods import REPO, completed, json_save, save, validate_boundary
from .rp1_14_prepare import CODONS, assert_fixed_mapping


def read_prepared(root: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict, dict]:
    status = json.loads((root / "completion.json").read_text())
    if not completed(root, status["fingerprint"]):
        raise ValueError("RP1-14 preparation is incomplete.")
    historical = pd.read_parquet(root / "reused_receptor_assignments.parquet")
    registry = pd.read_csv(root / "sample_registry.tsv", sep="\t")
    validate_boundary(historical, registry)
    annotation = historical[
        historical.productive_canonical & historical.rfu.notna()
    ].drop_duplicates("aminoAcid")
    annotation = annotation[["aminoAcid", "rfu", "max_cor", "pass_threshold"]]
    full, frames, coverages = {}, [], []
    for meta in registry.to_dict(orient="records"):
        path = root / "full_productive_counts" / f"{meta['sample_id']}.parquet"
        if not path.exists():
            if meta["status"] != "missing_raw_sample":
                raise ValueError("Prepared counts missing from an observed sample.")
            continue
        counts = pd.read_parquet(path)
        full[meta["sample_id"]] = counts.groupby("clone").read_count.sum()
        nts = counts.clone.str.split("|", n=1).str[0]
        counts["aminoAcid"] = [
            "".join(CODONS.get(s[i : i + 3], "X") for i in range(0, len(s), 3)) for s in nts
        ]
        mapped = counts.merge(annotation, on="aminoAcid", how="inner", validate="many_to_one")
        mapped["length_bin"] = mapped.aminoAcid.str.len() // 5
        mapped["vj"] = mapped.v_call + "|" + mapped.j_call
        mapped = mapped.assign(
            **{
                key: meta[key]
                for key in [
                    "sample_id",
                    "donor",
                    "compartment",
                    "visit",
                    "elapsed_years",
                    "dataset_id",
                ]
            }
        )
        frames.append(mapped)
        primary = mapped[mapped.pass_threshold]
        coverages.append(
            {
                **{
                    k: v
                    for k, v in meta.items()
                    if not k.startswith("source_") and k not in ["age_years", "authorized"]
                },
                "all_source_reads": meta["source_read_depth"],
                "fixed_nearest_reads": int(mapped.read_count.sum()),
                "fixed_nearest_receptors": mapped.clone.nunique(),
                "fixed_primary_reads": int(primary.read_count.sum()),
                "fixed_primary_receptors": primary.clone.nunique(),
                "fixed_primary_rfus": primary.rfu.nunique(),
                "fixed_primary_fraction_all_reads": float(
                    primary.read_count.sum() / meta["source_read_depth"]
                ),
                "fixed_primary_fraction_productive_reads": float(
                    primary.read_count.sum() / meta["productive_read_depth"]
                ),
                "fixed_primary_V_unresolved_read_fraction": float(
                    primary.loc[
                        primary.v_call.str.startswith("unresolved_or_family:"), "read_count"
                    ].sum()
                    / primary.read_count.sum()
                ),
                "fixed_primary_J_unresolved_read_fraction": float(
                    primary.loc[
                        primary.j_call.str.startswith("unresolved_or_family:"), "read_count"
                    ].sum()
                    / primary.read_count.sum()
                ),
            }
        )
    frame = pd.concat(frames, ignore_index=True)
    assert_fixed_mapping(frame)
    return frame, pd.DataFrame(coverages), full, status


def visit_counts(frame: pd.DataFrame, weighting: str) -> pd.Series:
    counts = frame.groupby("clone").read_count.sum()
    return counts if weighting == "reads" else counts.clip(upper=1)


def total_variation(a: pd.Series, b: pd.Series) -> float:
    if a.sum() <= 0 or b.sum() <= 0:
        return float("nan")
    return float((a / a.sum()).sub(b / b.sum(), fill_value=0).abs().sum() / 2)


def analysis_pairs(coverage: pd.DataFrame) -> pd.DataFrame:
    """Return observed chronological pairs; never impute missing samples."""
    records = []
    for (donor, compartment), group in coverage.groupby(["donor", "compartment"], sort=True):
        group = group.sort_values("visit")
        for a, b in itertools.combinations(group.to_dict(orient="records"), 2):
            records.append(
                {
                    "donor": donor,
                    "compartment": compartment,
                    "sample_before": a["sample_id"],
                    "sample_after": b["sample_id"],
                    "visit_before": a["visit"],
                    "visit_after": b["visit"],
                    "interval_years": b["elapsed_years"] - a["elapsed_years"],
                    "primary_pair": a["visit"] == group.visit.min()
                    and b["visit"] == group.visit.max(),
                }
            )
    return pd.DataFrame(records)


def persistence_rows(
    a: pd.DataFrame, b: pd.DataFrame, base: dict, detections: list[int]
) -> list[dict]:
    result = []
    ca = a.groupby(["rfu", "clone"]).read_count.sum()
    cb = b.groupby(["rfu", "clone"]).read_count.sum()
    rfu_a = set(a.rfu)
    rfu_b = set(b.rfu)
    for rfu in sorted(rfu_a | rfu_b):
        x = ca.loc[rfu] if rfu in rfu_a else pd.Series(dtype=float)
        y = cb.loc[rfu] if rfu in rfu_b else pd.Series(dtype=float)
        for detection in detections:
            result.append(
                {
                    **base,
                    "rfu": rfu,
                    "detection_reads": detection,
                    "before_reads": int(x.sum()),
                    "after_reads": int(y.sum()),
                    "observed_before": bool(x.sum() >= detection),
                    "observed_after": bool(y.sum() >= detection),
                    "clones_before": len(x),
                    "clones_after": len(y),
                    "shared_clones": len(x.index.intersection(y.index)),
                    "dominant_fraction_before": float(x.max() / x.sum()) if len(x) else np.nan,
                    "dominant_fraction_after": float(y.max() / y.sum()) if len(y) else np.nan,
                    "one_read_fraction_before": 1 / a.read_count.sum(),
                    "one_read_fraction_after": 1 / b.read_count.sum(),
                }
            )
    return result


def measurements(
    frame: pd.DataFrame, coverage: pd.DataFrame, full: dict, config: dict, out: Path
) -> pd.DataFrame:
    mapping = frame.drop_duplicates("clone").set_index("clone")
    samples = dict(tuple(frame.groupby("sample_id", sort=True)))
    coverage_index = coverage.set_index("sample_id")
    rows, support, persistence, full_rows, changes = [], [], [], [], []
    pairs = analysis_pairs(coverage)
    for base in pairs.to_dict(orient="records"):
        a, b = samples[base["sample_before"]], samples[base["sample_after"]]
        full_rows.append(
            {
                **base,
                "productive_clone_TV_reads": total_variation(
                    full[base["sample_before"]], full[base["sample_after"]]
                ),
                "productive_clone_TV_unique": total_variation(
                    full[base["sample_before"]].clip(upper=1),
                    full[base["sample_after"]].clip(upper=1),
                ),
            }
        )
        for policy in ("threshold", "nearest"):
            x, y = (d[d.pass_threshold] if policy == "threshold" else d for d in (a, b))
            ok = min(x.clone.nunique(), y.clone.nunique()) >= config["min_receptors_per_visit"]
            support.append(
                {
                    **base,
                    "policy": policy,
                    "receptors_before": x.clone.nunique(),
                    "receptors_after": y.clone.nunique(),
                    "status": "analyzed" if ok else "insufficient_observed_receptors",
                }
            )
            if not ok:
                continue
            cov = {
                f"retained_fraction_productive_reads_{side}": float(
                    d.read_count.sum()
                    / coverage_index.loc[base[f"sample_{side}"], "productive_read_depth"]
                )
                for side, d in [("before", x), ("after", y)]
            }
            for weight in ("reads", "unique_clone"):
                ca, cb = visit_counts(x, weight), visit_counts(y, weight)
                for feature, grouping in [("rfu", "RFU"), ("v_call", "TRBV"), ("vj", "TRBV_TRBJ")]:
                    change = multiscale_repertoire_change(ca, cb, mapping[feature])
                    record = {
                        **base,
                        **cov,
                        "policy": policy,
                        "weighting": weight,
                        "grouping": grouping,
                        **change.summary,
                    }
                    rows.append(record)
                    if grouping == "RFU" and policy == "threshold" and weight == "reads":
                        changes.append(change.group_changes.reset_index().assign(**base))
            if policy == "threshold":
                persistence.extend(persistence_rows(x, y, base, config["rfu_detection_reads"]))
    table = pd.DataFrame(rows)
    save(table, out, "rp1_14_multiscale_pairs.tsv")
    save(pd.DataFrame(support), out, "rp1_14_pair_support.tsv")
    save(pd.DataFrame(persistence), out, "rp1_14_rfu_persistence.tsv")
    save(pd.DataFrame(full_rows), out, "rp1_14_full_productive_clone_distances.tsv")
    save(pd.concat(changes, ignore_index=True), out, "rp1_14_rfu_changes.tsv")
    primary = table[
        table.primary_pair
        & table.policy.eq("threshold")
        & table.weighting.eq("reads")
        & table.grouping.eq("RFU")
    ]
    paired = primary.pivot(
        index=["donor", "visit_before", "visit_after"],
        columns="compartment",
        values=["d_clone", "d_group", "aggregation_cancellation"],
    )
    paired.columns = [f"{metric}_{comp}" for metric, comp in paired.columns]
    paired = paired.dropna().reset_index()
    paired["RFU_TV_CD8_minus_CD4"] = paired.d_group_CD8 - paired.d_group_CD4
    save(paired, out, "rp1_14_paired_compartment_comparison.tsv")
    return table


def similarities(frame: pd.DataFrame, out: Path) -> None:
    tables, summaries = [], []
    selected = frame[frame.pass_threshold]
    for compartment, part in selected.groupby("compartment"):
        matrix = part.pivot_table(
            index="sample_id", columns="rfu", values="read_count", aggfunc="sum", fill_value=0
        )
        meta = part[["sample_id", "donor", "visit"]].drop_duplicates().set_index("sample_id")
        pairs = longitudinal_similarity(
            matrix, metadata=meta, donor_key="donor", time_key="visit", metric="cosine"
        )
        tables.append(pairs.assign(compartment=compartment))
        for donor in sorted(meta.donor.unique()):
            for within in (True, False):
                rows = pairs[
                    pairs.same_donor.eq(within)
                    & (pairs.donor_a.eq(donor) | pairs.donor_b.eq(donor))
                ]
                summaries.append(
                    {
                        "donor": donor,
                        "compartment": compartment,
                        "relation": "within" if within else "between_mean",
                        "mean_cosine": float(rows.value.mean()),
                        "n_dependent_pairs": len(rows),
                    }
                )
    save(pd.concat(tables, ignore_index=True), out, "rp1_14_similarity_pairs.tsv")
    save(pd.DataFrame(summaries), out, "rp1_14_similarity_donor_summary.tsv")


def empirical_reads(counts: pd.Series, n: int, rng: np.random.Generator) -> pd.Series:
    """Hypergeometric draw of actual observed reads, never multinomial/templates."""
    values = counts.to_numpy()
    if (
        not np.isfinite(values).all()
        or (values < 0).any()
        or not np.equal(values, np.floor(values)).all()
    ):
        raise ValueError("Empirical read subsampling requires actual integer counts.")
    if n > counts.sum() or n < 0:
        raise ValueError("Read subsampling exceeds observed depth.")
    draw = pd.Series(
        rng.multivariate_hypergeometric(values.astype(np.int64), n), index=counts.index
    )
    return draw[draw > 0]


def sensitivities(frame: pd.DataFrame, coverage: pd.DataFrame, config: dict, out: Path) -> None:
    selected = frame[frame.pass_threshold]
    mapping = selected.drop_duplicates("clone").set_index("clone")
    strata = mapping[["v_call", "j_call", "length_bin"]]
    exchangeable = float(
        strata.assign(rfu=mapping.rfu).groupby(list(strata)).rfu.transform("nunique").gt(1).mean()
    )
    pairs = analysis_pairs(coverage)
    pairs = pairs[pairs.primary_pair]
    sample_counts = {s: visit_counts(g, "reads") for s, g in selected.groupby("sample_id")}
    controls, sampling, dominance = [], [], []
    for replicate in range(config["random_group_replicates"]):
        grouping = permute_fixed_groups(
            mapping.rfu, strata=strata, random_state=config["random_seed"] + replicate
        )
        changed = float(grouping.ne(mapping.rfu).mean())
        for base in pairs.to_dict(orient="records"):
            a, b = sample_counts[base["sample_before"]], sample_counts[base["sample_after"]]
            change = multiscale_repertoire_change(a, b, grouping)
            controls.append(
                {
                    **base,
                    "replicate": replicate,
                    "seed": config["random_seed"] + replicate,
                    "exchangeable_receptor_fraction": exchangeable,
                    "changed_label_fraction": changed,
                    **change.summary,
                }
            )
    for i, base in enumerate(pairs.to_dict(orient="records")):
        a, b = sample_counts[base["sample_before"]], sample_counts[base["sample_after"]]
        n = int(min(config["depth_cap_reads"], a.sum(), b.sum()))
        for rep in range(config["depth_replicates"]):
            seed = config["random_seed"] + 10000 * (i + 1) + rep
            rng = np.random.default_rng(seed)
            x, y = empirical_reads(a, n, rng), empirical_reads(b, n, rng)
            change = multiscale_repertoire_change(x, y, mapping.rfu)
            sampling.append(
                {
                    **base,
                    "replicate": rep,
                    "seed": seed,
                    "reads_per_visit": n,
                    "scheme": "observed_reads_without_replacement_conditional_on_reuse_universe",
                    **change.summary,
                }
            )
        removed = {a.idxmax(), b.idxmax()}
        x, y = a.drop(list(removed), errors="ignore"), b.drop(list(removed), errors="ignore")
        change = multiscale_repertoire_change(x, y, mapping.rfu)
        dominance.append(
            {
                **base,
                "removed_receptors": len(removed),
                "retained_read_fraction_before": float(x.sum() / a.sum()),
                "retained_read_fraction_after": float(y.sum() / b.sum()),
                **change.summary,
            }
        )
    save(pd.DataFrame(controls), out, "rp1_14_fixed_group_controls.tsv")
    save(pd.DataFrame(sampling), out, "rp1_14_empirical_read_sensitivity.tsv")
    save(pd.DataFrame(dominance), out, "rp1_14_dominant_clone_sensitivity.tsv")


def summarize(out: Path, prepared: Path) -> dict:
    counts = pd.read_csv(out / "rp1_14_sample_coverage.tsv", sep="\t")
    table = pd.read_csv(out / "rp1_14_multiscale_pairs.tsv", sep="\t")
    primary = table[table.primary_pair & table.policy.eq("threshold") & table.weighting.eq("reads")]
    p = pd.read_csv(out / "rp1_14_rfu_persistence.tsv", sep="\t")
    p = p[p.primary_pair & p.detection_reads.eq(5)]
    pr = []
    for (donor, comp), x in p.groupby(["donor", "compartment"]):
        persistent = x[x.observed_before & x.observed_after]
        pr.append(
            {
                "donor": donor,
                "compartment": comp,
                "persistent_rfus": len(persistent),
                "union_detected_rfus": int((x.observed_before | x.observed_after).sum()),
                "persistent_fraction_union": len(persistent)
                / (x.observed_before | x.observed_after).sum(),
                "persistent_without_shared_clone": int(persistent.shared_clones.eq(0).sum()),
                "fraction_persistent_without_shared_clone": float(
                    persistent.shared_clones.eq(0).mean()
                ),
                "median_dominant_fraction_before": float(
                    persistent.dominant_fraction_before.median()
                ),
                "median_dominant_fraction_after": float(
                    persistent.dominant_fraction_after.median()
                ),
            }
        )
    save(pd.DataFrame(pr), out, "rp1_14_persistence_donor_summary.tsv")
    summary = {
        "samples": len(counts),
        "donors": counts.donor.nunique(),
        "source_rows": int(
            pd.read_csv(prepared / "sample_registry.tsv", sep="\t").source_rows.sum()
        ),
        "source_reads": int(counts.all_source_reads.sum()),
        "productive_reads": int(counts.productive_read_depth.sum()),
        "primary_mapped_reads": int(counts.fixed_primary_reads.sum()),
        "primary_sample_receptor_observations": int(counts.fixed_primary_receptors.sum()),
        "primary_coverage_all_reads_min_max": [
            float(counts.fixed_primary_fraction_all_reads.min()),
            float(counts.fixed_primary_fraction_all_reads.max()),
        ],
        "primary_coverage_productive_reads_min_max": [
            float(counts.fixed_primary_fraction_productive_reads.min()),
            float(counts.fixed_primary_fraction_productive_reads.max()),
        ],
        "interval_years_min_max": [
            float(primary.interval_years.min()),
            float(primary.interval_years.max()),
        ],
        "distance_rows": len(table),
        "tv_inequality_failures": int((table.d_group > table.d_clone + 1e-12).sum()),
        "primary_by_compartment_grouping": [],
    }
    for (comp, group), x in primary.groupby(["compartment", "grouping"]):
        summary["primary_by_compartment_grouping"].append(
            {
                "compartment": comp,
                "grouping": group,
                "donors": x.donor.nunique(),
                "clone_tv_median": float(x.d_clone.median()),
                "group_tv_median": float(x.d_group.median()),
                "cancellation_median": float(x.aggregation_cancellation.median()),
            }
        )
    return summary


def figures(out: Path, fingerprint: str, command: str) -> None:
    import matplotlib.pyplot as plt

    set_style()
    coverage = pd.read_csv(out / "rp1_14_sample_coverage.tsv", sep="\t")
    table = pd.read_csv(out / "rp1_14_multiscale_pairs.tsv", sep="\t")
    primary = table[
        table.primary_pair
        & table.policy.eq("threshold")
        & table.weighting.eq("reads")
        & table.grouping.eq("RFU")
    ]
    similarity = pd.read_csv(out / "rp1_14_similarity_donor_summary.tsv", sep="\t")
    persistence = pd.read_csv(out / "rp1_14_persistence_donor_summary.tsv", sep="\t")
    controls = pd.read_csv(out / "rp1_14_fixed_group_controls.tsv", sep="\t")
    donors = sorted(coverage.donor.unique())
    colors = {"CD4": "#0072B2", "CD8": "#D55E00"}
    fig, axs = plt.subplots(2, 3, figsize=(10.5, 6.6), constrained_layout=True)
    ax = axs[0, 0]
    for comp, offset in [("CD4", -0.12), ("CD8", 0.12)]:
        for donor, part in coverage[coverage.compartment.eq(comp)].groupby("donor"):
            y = donors.index(donor) + offset
            ax.plot(part.elapsed_years, [y] * len(part), color=colors[comp], alpha=0.45)
            ax.scatter(
                part.elapsed_years,
                [y] * len(part),
                s=20 + 100 * part.fixed_primary_fraction_productive_reads,
                facecolors="none" if comp == "CD8" else colors[comp],
                edgecolors=colors[comp],
                label=comp if donor == donors[0] else None,
            )
    ax.set(
        yticks=range(len(donors)),
        yticklabels=donors,
        xlabel="Years from first collection",
        title="Design and reuse coverage",
    )
    ax.legend(
        fontsize=6,
        title="Size: qualified fraction\nof productive reads",
        title_fontsize=6,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.23),
        ncol=2,
    )
    ax = axs[0, 1]
    for comp, offset in [("CD4", -0.08), ("CD8", 0.08)]:
        for _, part in similarity[similarity.compartment.eq(comp)].groupby("donor"):
            values = part.set_index("relation").mean_cosine
            ax.plot(
                [0 + offset, 1 + offset],
                [values["within"], values["between_mean"]],
                "o-",
                color=colors[comp],
                alpha=0.7,
            )
    ax.set(
        xticks=[0, 1],
        xticklabels=["Within donor", "Between donors"],
        ylabel="RFU cosine similarity",
        title="Donor-specific RFU profiles",
        ylim=(0, 1),
    )
    ax = axs[0, 2]
    for comp, part in primary.groupby("compartment"):
        ax.scatter(part.d_clone, part.d_group, color=colors[comp], label=comp)
    ax.plot([0, 1], [0, 1], color=".6", lw=0.7, ls="--")
    ax.set(
        xlabel="Nucleotide/V/J receptor TV",
        ylabel="RFU TV",
        title="Change hidden by aggregation",
        xlim=(0, 1),
        ylim=(0, 1),
    )
    ax.legend()
    ax = axs[1, 0]
    for _, part in primary.groupby("donor"):
        x = part.set_index("compartment").d_group
        ax.plot([0, 1], [x["CD4"], x["CD8"]], "o-", color=".4", alpha=0.7)
    ax.set(
        xticks=[0, 1],
        xticklabels=["CD4", "CD8"],
        ylabel="RFU TV",
        title="Paired long-term remodeling",
        ylim=(0, 1),
    )
    ax = axs[1, 1]
    for comp, part in persistence.groupby("compartment"):
        ax.scatter(
            part.persistent_fraction_union,
            part.fraction_persistent_without_shared_clone,
            color=colors[comp],
            label=comp,
        )
    ax.set(
        xlabel="Persistent / union detected RFUs",
        ylabel="Fraction with no shared observed receptor",
        title="Persistent groups, changing receptors",
        xlim=(0, 1),
        ylim=(0, 1),
    )
    ax.text(0.03, 0.03, "Detection ≥5 reads per visit", fontsize=6, transform=ax.transAxes)
    ax = axs[1, 2]
    for comp, offset in [("CD4", -0.13), ("CD8", 0.13)]:
        for i, donor in enumerate(donors):
            x = controls[
                controls.donor.eq(donor) & controls.compartment.eq(comp)
            ].aggregation_cancellation
            obs = primary[
                primary.donor.eq(donor) & primary.compartment.eq(comp)
            ].aggregation_cancellation.iloc[0]
            lo, med, hi = x.quantile([0.025, 0.5, 0.975])
            ax.plot([i + offset] * 2, [lo, hi], color=colors[comp], alpha=0.55)
            ax.scatter(i + offset, med, facecolor="white", edgecolor=colors[comp], s=20)
            ax.scatter(i + offset, obs, color=colors[comp], marker="D", s=15)
    ax.set(
        xticks=range(len(donors)),
        xticklabels=donors,
        ylabel="Aggregation cancellation",
        title="Fixed feature-matched controls",
        ylim=(0, 0.75),
    )
    ax.text(
        0.03,
        0.97,
        "◆ RFU   ○ matched-map median\nLines: 30-map percentile range",
        va="top",
        fontsize=6,
        transform=ax.transAxes,
    )
    for letter, ax in zip("ABCDEF", axs.flat, strict=True):
        clean_axis(ax)
        panel_label(ax, letter)
    save_figure(fig, out, "rp1_14_longitudinal")
    plt.close(fig)
    definitions = [
        (
            "A",
            "rp1_14_sample_coverage.tsv",
            "six donors; 36 samples; productive-read denominator",
            "Historical top-row reuse dictionary; read coverage is not template/cell coverage",
        ),
        (
            "B",
            "rp1_14_similarity_donor_summary.tsv",
            "one within/between average per donor/compartment",
            "All observed visit pairs; correlated averages, no pair-level hypothesis test",
        ),
        (
            "C",
            "rp1_14_multiscale_pairs.tsv",
            "six donors per compartment; earliest/latest; threshold reads",
            "Same fixed reuse universe for receptor and group TV; excluded mass reported",
        ),
        (
            "D",
            "rp1_14_paired_compartment_comparison.tsv",
            "six CD4/CD8 paired donors",
            "Read-weighted relative change; small descriptive cohort",
        ),
        (
            "E",
            "rp1_14_persistence_donor_summary.tsv",
            "six donors per compartment; RFUs detected at >=5 reads at both visits",
            "Nondetection is not absence; persistence is not antigen-functional stability",
        ),
        (
            "F",
            "rp1_14_fixed_group_controls.tsv",
            "six donors per compartment; 30 repeated fixed group maps",
            "Descriptive granularity/feature control; not calibrated null pvalues",
        ),
    ]
    save(
        pd.DataFrame(
            [
                {
                    "figure": "rp1_14_longitudinal.pdf",
                    "panel": p,
                    "source_table": s,
                    "source_sha256": file_sha256(out / s),
                    "denominator": d,
                    "limits": lim,
                    "fingerprint": fingerprint,
                    "command": command,
                    "status": "complete",
                }
                for p, s, d, lim in definitions
            ]
        ),
        out,
        "figure_source_index.tsv",
    )


def validate_measurement_stage(out: Path, identity: dict) -> bool:
    path = out / "measurement_stage.json"
    if not path.exists():
        return False
    stage = json.loads(path.read_text())
    if stage.get("status") != "complete" or stage.get("identity") != identity:
        raise ValueError("Incomplete measurement stage or changed scientific inputs/settings.")
    for name, checksum in stage["outputs"].items():
        if not (out / name).is_file() or file_sha256(out / name) != checksum:
            raise ValueError("Measurement-stage table checksum failed.")
    return True


def run(args: argparse.Namespace) -> None:
    out = args.out.resolve()
    if out.is_relative_to(REPO):
        raise ValueError("Biological results must remain outside Git.")
    config = json.loads(args.config.read_text())
    prepared = json.loads((args.prepared / "completion.json").read_text())
    if not completed(args.prepared, prepared["fingerprint"]):
        raise ValueError("Preparation incomplete.")
    identity = {
        "preparation_manifest_sha256": file_sha256(args.prepared / "completion.json"),
        "config_sha256": file_sha256(args.config),
        "code_sha256": {
            str(p.relative_to(REPO)): file_sha256(p)
            for p in [
                Path(__file__),
                REPO / "src/scrfu/multiscale.py",
                REPO / "src/scrfu/longitudinal.py",
                REPO / "manuscript/scripts/_figure_common.py",
            ]
        },
        "python": sys.executable,
        "scrfu_import": scrfu.__file__,
        "versions": {
            k: importlib.metadata.version(k) for k in ["numpy", "pandas", "matplotlib", "scrfu"]
        },
    }
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    if not args.refresh_exports and completed(out, fingerprint):
        print("Completed RP1-14 longitudinal analysis verified; no measurements rerun.")
        return
    measurement_identity = json.loads(json.dumps(identity))
    measurement_identity["code_sha256"].pop(str(Path(__file__).relative_to(REPO)))
    measurement_identity["measurement_functions"] = {
        name: hashlib.sha256(inspect.getsource(globals()[name]).encode()).hexdigest()
        for name in [
            "read_prepared",
            "visit_counts",
            "total_variation",
            "analysis_pairs",
            "persistence_rows",
            "measurements",
            "similarities",
            "empirical_reads",
            "sensitivities",
        ]
    }
    stage_path = out / "measurement_stage.json"
    reused = validate_measurement_stage(out, measurement_identity)
    if (out / "completion.json").exists() and not reused:
        previous = json.loads((out / "completion.json").read_text())
        if previous.get("fingerprint") != fingerprint:
            raise ValueError(
                "Partial analysis settings changed: use a new versioned output directory."
            )
    out.mkdir(parents=True, exist_ok=True)
    json_save({"status": "running", "fingerprint": fingerprint}, out / "completion.json")
    frame, coverage, full, _ = read_prepared(args.prepared)
    if not reused:
        save(coverage, out, "rp1_14_sample_coverage.tsv")
        print(
            f"Fixed reuse dictionary mapped {len(frame)} sample/receptor observations; {frame.donor.nunique()} donors",
            flush=True,
        )
        measurements(frame, coverage, full, config, out)
        print("Paired distances and persistence complete", flush=True)
        similarities(frame, out)
        sensitivities(frame, coverage, config, out)
        json_save(
            {
                "status": "complete",
                "identity": measurement_identity,
                "outputs": {p.name: file_sha256(p) for p in out.glob("*.tsv")},
            },
            stage_path,
        )
    else:
        print("Verified completed measurement stage reused; resuming export only", flush=True)
    summary = summarize(out, args.prepared)
    summary.update(
        {
            "primary_unique_receptors": frame.loc[frame.pass_threshold, "clone"].nunique(),
            "primary_unique_rfus": frame.loc[frame.pass_threshold, "rfu"].nunique(),
        }
    )
    json_save(summary, out / "evidence_counts.json")
    figures(out, fingerprint, " ".join(sys.argv))
    json_save(
        {
            **identity,
            "fingerprint": fingerprint,
            "config": config,
            "command": sys.argv,
            "git_sha": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
            ).strip(),
            "source_scope": "published_authorized_RP1_only",
            "assignment": "historical_vectors_reused_after_verified_alignment_repair",
            "fixed_universe": "reused_AA_RFUs_applied_to_full_productive_counts_each_visit",
            "historical_reference_identity": "bounded_parity_only; original execution provenance unavailable",
        },
        out / "provenance.json",
    )
    outputs = {
        p.name: file_sha256(p)
        for p in sorted(out.iterdir())
        if p.is_file() and p.name != "completion.json"
    }
    json_save(
        {"status": "complete", "fingerprint": fingerprint, "outputs": outputs},
        out / "completion.json",
    )
    print(json.dumps(summary, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--refresh-exports",
        action="store_true",
        help="Refresh exports only after measurement input/function hashes pass; reuse completed measurements.",
    )
    parser.add_argument("--config", type=Path, default=REPO / "manuscript/config/rp1_14_v1.json")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
