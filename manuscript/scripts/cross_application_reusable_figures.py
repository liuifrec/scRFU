"""Data-only exports from completed scRFU applications; no biological runners."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import shlex
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

from .rp1_14_reusable_figures import REPO, json_save, require, safe_child, sha256

CONFIG = REPO / "manuscript/config/cross_application_reusable_figures_v1.json"
CAPTIONS = REPO / "manuscript/figures/cross_application_reusable_captions.json"
SCHEMA = REPO / "manuscript/figures/cross_application_table_schema.json"
SOURCES = {
    "a_visits": ("development", "gse190905_paired_manifest.tsv"),
    "a_distances": ("development", "gse190905_multiscale_pairs.tsv"),
    "a_support": ("development", "gse190905_pair_support.tsv"),
    "a_depth": ("development", "gse190905_empirical_depth_sensitivity.tsv"),
    "a_dominance": ("development", "gse190905_dominant_clone_sensitivity.tsv"),
    "aliases": ("development", "gse190905_source_donor_crosswalk.tsv"),
    "dev_evidence": ("development", "evidence_counts.json"),
    "b_visits": ("external", "visit_support.tsv"),
    "b_pairs": ("external", "pair_support.tsv"),
    "b_distances": ("external", "multiscale_pairs.tsv"),
    "b_persistence": ("external", "persistence_summary.tsv"),
    "b_evidence": ("external", "evidence_counts.json"),
    "cross": ("finish", "cross_application_summary.tsv"),
    "rp_evidence": ("rp", "evidence_counts.json"),
    "reg_evidence": ("regulatory", "evidence_counts.json"),
    "reg_candidates": ("regulatory", "high_confidence_regulatory_candidates.tsv"),
    "reg_targets": ("regulatory", "high_confidence_independent_qtl_chains.tsv"),
    "ordinary_deltas": ("ordinary", "paired_donor_deltas.tsv"),
    "ordinary_metrics": ("ordinary", "donor_held_out_metrics.tsv"),
    "purged_deltas": ("purged", "paired_donor_deltas.tsv"),
    "purged_metrics": ("purged", "donor_held_out_metrics.tsv"),
    "w_coverage": ("development", "wells_donor_tissue_coverage.tsv"),
    "w_support": ("development", "wells_rfu_support.tsv"),
    "w_context": ("development", "wells_rfu_donor_tissue_celltype.tsv"),
}
INPUTS_BY_FIGURE = {
    "A": ["a_visits", "a_distances", "a_support", "a_depth", "a_dominance", "dev_evidence"],
    "B": ["b_visits", "b_pairs", "b_distances", "b_persistence", "b_evidence"],
    "C": ["cross", "rp_evidence", "a_distances", "b_distances"],
    "D": ["reg_evidence", "reg_candidates", "reg_targets"],
    "E": ["ordinary_deltas", "ordinary_metrics", "purged_deltas", "purged_metrics", "aliases"],
    "F": ["w_coverage", "w_support", "w_context", "dev_evidence"],
}
RULES = {
    "A": "threshold policy, cell weights, RFU grouping; subsets all_T/CD4/CD8. Six paired donors, >=100 qualified cells/visit. Sensitivity: unique_clone, 50 frozen observed-cell draws at min(500,pair depth), dominant-receptor union removal. No resampling performed.",
    "B": "threshold policy, cell weights, RFU grouping. Preserve all 18 registered visits and 18 interval combinations; only 11 visits/8 pairs supported. Primary persistence at one cell/endpoint; frozen five-cell sensitivity retained, undefined fractions stay missing.",
    "C": "Copy the seven frozen cross-application summary rows; compare with upstream medians without new endpoints. Keep read/cell/GEX-matched-cell units, intervals, donor counts and distinct persistence rules. No pooling.",
    "D": "Use frozen source conditional-overlap counts; tier A_same_variant_conditional_layers for six-variant matrix. Deduplicate molecular target/rank records across RFU joins. Preserve direct_TCR_gene/non_TCR_gene classes. No QTL search, causal arrows or colocalization claim.",
    "E": "Clone-weighted primary log loss from six saved held-out folds in each design. Extended minus baseline; positive=worse. Equal-donor means, not cell pooling. Retain cell-weighted rows in source table. Frozen patient_tcr-to-donor aliases; no fit or probability rescoring.",
    "F": "Aggregate frozen qualified RFU context by source tissue/cell_type. Keep all 24 atlas donors and 10 tissues; missing atlas combinations are distinct from observed zero-TRB samples. Frozen support >=50 cells across >=4 donors; no assignment or specificity inference.",
}
DIST = ["d_clone", "d_group", "aggregation_cancellation"]
PAIR_B = ["donor", "compartment", "interval"]
INTERVALS = ["pre_to_last_radiation", "last_radiation_to_6weeks", "pre_to_6weeks"]


def verify_sources(workspaces: dict, config: dict) -> tuple[dict, dict, dict]:
    roots, verified = {}, {}
    for key, spec in config["upstream"].items():
        root = workspaces[spec["workspace"]] / spec["directory"]
        path = root / "completion.json"
        require(path.is_file() and sha256(path) == spec["sha256"], f"Changed completion: {path}")
        completed = json.loads(path.read_text())
        require(completed["status"] == "complete", f"Incomplete source: {root}")
        for name, checksum in completed["outputs"].items():
            output = safe_child(root, name)
            require(
                output.is_file() and sha256(output) == checksum, f"Upstream hash failed: {output}"
            )
        roots[key] = root
        verified[key] = {
            "path": str(path),
            "sha256": spec["sha256"],
            "verified_outputs": len(completed["outputs"]),
        }
    data, inputs = {}, {}
    for key, (stage, name) in SOURCES.items():
        path = roots[stage] / name
        data[key] = (
            json.loads(path.read_text()) if path.suffix == ".json" else pd.read_csv(path, sep="\t")
        )
        inputs[key] = {"path": str(path), "sha256": sha256(path), "upstream": stage}
    return data, inputs, verified


def primary(frame: pd.DataFrame, weights: str = "cell") -> pd.DataFrame:
    out = frame[
        frame.policy.eq("threshold") & frame.weighting.eq(weights) & frame.grouping.eq("RFU")
    ].copy()
    require(out.status.eq("valid").all(), "Invalid frozen primary measurement.")
    require(
        np.allclose(out.d_clone - out.d_group, out.aggregation_cancellation),
        "Cancellation identity failed.",
    )
    return out


def sensitivity_table(primary_a, distances, depth, dominance) -> pd.DataFrame:
    """Use within-donor technical medians, retaining frozen ranges and depths."""
    rows = []
    for condition, frame in [
        ("cell", primary_a),
        ("unique", primary(distances, "unique_clone").query("subset == 'all_T'")),
    ]:
        require(len(frame) == frame.donor.nunique() == 6, "Incomplete sensitivity donors.")
        for r in frame.itertuples():
            rows.append(
                dict(
                    donor=r.donor,
                    condition=condition,
                    rfu_tv=r.d_group,
                    low=r.d_group,
                    high=r.d_group,
                    replicates=1,
                    cells_per_visit=np.nan,
                    retained_fraction_before=1.0,
                    retained_fraction_after=1.0,
                )
            )
    require(
        len(depth) == 300 and not depth.duplicated(["donor", "replicate"]).any(),
        "Depth replicate coverage differs.",
    )
    for donor, part in depth.groupby("donor"):
        require(
            set(part.replicate) == set(range(50)) and part.cells_per_visit.nunique() == 1,
            "Incomplete technical draws.",
        )
        lo, med, hi = part.d_group.quantile([0.025, 0.5, 0.975], interpolation="linear")
        rows.append(
            dict(
                donor=donor,
                condition="depth",
                rfu_tv=med,
                low=lo,
                high=hi,
                replicates=50,
                cells_per_visit=part.cells_per_visit.iloc[0],
                retained_fraction_before=np.nan,
                retained_fraction_after=np.nan,
            )
        )
    require(
        sorted(depth.groupby("donor").cells_per_visit.first()) == [404, 500, 500, 500, 500, 500],
        "Frozen common depths differ.",
    )
    for r in dominance.itertuples():
        rows.append(
            dict(
                donor=r.donor,
                condition="dominant_removed",
                rfu_tv=r.d_group,
                low=r.d_group,
                high=r.d_group,
                replicates=1,
                cells_per_visit=np.nan,
                retained_fraction_before=r.retained_cell_fraction_before,
                retained_fraction_after=r.retained_cell_fraction_after,
            )
        )
    out = pd.DataFrame(rows)
    require(
        len(out) == 24 and not out.duplicated(["donor", "condition"]).any(),
        "Sensitivity pairing differs.",
    )
    return out.sort_values(["donor", "condition"]).reset_index(drop=True)


def external_pairs(support: pd.DataFrame, distance: pd.DataFrame) -> pd.DataFrame:
    require(
        not support.duplicated(PAIR_B).any() and not distance.duplicated(PAIR_B).any(),
        "Duplicate external pair.",
    )
    result = support.merge(distance[[*PAIR_B, *DIST]], on=PAIR_B, how="left", validate="one_to_one")
    require(
        result.status.eq("analyzed").eq(result.d_group.notna()).all(),
        "Unavailable interval acquired an estimate.",
    )
    return result


def regulatory_matrix(
    candidates: pd.DataFrame, targets: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    tier = "A_same_variant_conditional_layers"
    selected = candidates[candidates.evidence_tier.eq(tier)].copy()
    require(
        len(selected) == len(selected.drop_duplicates(["variant", "rfu_label"])) == 62,
        "Variant/RFU pair count differs.",
    )
    require(
        selected.variant.nunique() == 6 and selected.rfu_label.nunique() == 36,
        "Same-variant overlap differs.",
    )
    records = (
        targets[targets.evidence_tier.eq(tier)].drop(columns=["rfu_label"]).drop_duplicates().copy()
    )
    require(
        not records.duplicated(["variant", "layer", "target", "conditional_rank"]).any(),
        "Conflicting molecular records.",
    )
    require(
        records.layer.value_counts().to_dict() == {"caqtl": 10, "eqtl": 6},
        "Molecular records double-counted across RFUs.",
    )
    rows = []
    for variant, part in records.groupby("variant", sort=True):
        eq = part[part.layer.eq("eqtl")]
        ca = part[part.layer.eq("caqtl")]
        require(
            len(eq) == 1 and eq.target_class.isin(["direct_TCR_gene", "non_TCR_gene"]).all(),
            "Source target classes changed.",
        )
        rows.append(
            dict(
                variant=variant,
                eqtl_target=eq.target.iloc[0],
                target_class=eq.target_class.iloc[0],
                eqtl_records=len(eq),
                caqtl_records=len(ca),
                caqtl_targets=";".join(sorted(ca.target.unique())),
                rfu_count=selected[selected.variant.eq(variant)].rfu_label.nunique(),
            )
        )
    return pd.DataFrame(rows), records.sort_values(
        ["variant", "layer", "target", "conditional_rank"]
    ).reset_index(drop=True)


def prediction_tables(data: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    alias = data["aliases"].set_index("patient_tcr").donor
    rows = []
    for design in ["ordinary", "purged"]:
        part = data[f"{design}_deltas"].copy()
        metrics = data[f"{design}_metrics"]
        require(
            not part.duplicated(["held_out_donor", "weighting"]).any(), "Duplicate prediction fold."
        )
        require(
            np.allclose(part.extended_log_loss - part.baseline_log_loss, part.delta_log_loss),
            "Prediction delta sign differs.",
        )
        pivot = metrics.pivot(
            index=["held_out_donor", "weighting"], columns="model", values="log_loss"
        )
        indexed = part.set_index(["held_out_donor", "weighting"]).sort_index()
        for col, model in [
            ("baseline_log_loss", "TRBV_TRBJ_length"),
            ("extended_log_loss", "TRBV_TRBJ_length_RFU"),
        ]:
            require(
                np.allclose(indexed[col], pivot[model].sort_index()),
                "Saved metrics and deltas disagree.",
            )
        part["donor"] = part.held_out_donor.map(alias)
        require(
            part.donor.notna().all() and part.donor.nunique() == 6, "Frozen donor alias missing."
        )
        part = part.drop(columns="held_out_donor").assign(design=design)
        rows.append(part)
    all_rows = (
        pd.concat(rows, ignore_index=True)
        .sort_values(["design", "weighting", "donor"])
        .reset_index(drop=True)
    )
    selected = all_rows[all_rows.weighting.eq("clone")]
    require(
        len(selected) == 12 and selected.delta_log_loss.gt(0).all(),
        "Frozen primary prediction direction differs.",
    )
    summary = selected.groupby("design", as_index=False).agg(
        donors=("donor", "nunique"),
        baseline_log_loss=("baseline_log_loss", "mean"),
        extended_log_loss=("extended_log_loss", "mean"),
        delta_log_loss=("delta_log_loss", "mean"),
    )
    for design, expected in [
        ("ordinary", [1.850996, 1.855212, 0.004216]),
        ("purged", [1.852679, 1.856996, 0.004317]),
    ]:
        row = summary[summary.design.eq(design)].iloc[0]
        require(
            np.allclose(
                row[["baseline_log_loss", "extended_log_loss", "delta_log_loss"]].to_numpy(
                    dtype=float
                ),
                expected,
                rtol=0,
                atol=5e-7,
            ),
            "Prediction means disagree with frozen report.",
        )
    return all_rows, summary


def wells_grid(coverage: pd.DataFrame) -> pd.DataFrame:
    require(not coverage.duplicated(["donor", "tissue"]).any(), "Duplicate Wells donor/tissue.")
    index = pd.MultiIndex.from_product(
        [sorted(coverage.donor.unique()), sorted(coverage.tissue.unique())],
        names=["donor", "tissue"],
    )
    grid = coverage.set_index(["donor", "tissue"]).reindex(index).reset_index()
    grid["coverage_state"] = np.select(
        [grid.atlas_cells.isna(), grid.primary_trb_cells.eq(0)],
        ["no_atlas_sample", "atlas_zero_primary_TRB"],
        default="primary_TRB_observed",
    )
    require(
        grid.loc[grid.primary_trb_cells.eq(0), "threshold_fraction_of_primary_trb"].isna().all(),
        "Zero-TRB coverage must remain undefined.",
    )
    require(
        np.allclose(
            grid.threshold_cells / grid.primary_trb_cells.replace(0, np.nan),
            grid.threshold_fraction_of_primary_trb,
            equal_nan=True,
        ),
        "Wells coverage denominator differs.",
    )
    return grid


def build_tables(data: dict) -> dict:
    visits = data["a_visits"].copy()
    require(len(visits) == 12 and visits.donor.nunique() == 6, "Cached GSE190905 design differs.")
    require(
        visits.tcr_cells.sum() == 27655 and visits.threshold_cells.sum() == 21657,
        "GSE190905 cell counts differ.",
    )
    require(
        np.allclose(visits.threshold_cells / visits.tcr_cells, visits.threshold_coverage_of_tcr),
        "GSE190905 coverage denominator differs.",
    )
    visits = (
        visits.drop(columns=["authorized", "library"])
        .sort_values(["donor", "time"])
        .reset_index(drop=True)
    )
    a = primary(data["a_distances"])
    a = a[a.subset.isin(["all_T", "CD4", "CD8"])].copy()
    require(
        a.groupby("subset").donor.nunique().eq(6).all() and len(a) == 18,
        "GSE190905 subset support differs.",
    )
    a = (
        a[
            [
                "donor",
                "subset",
                "treatment",
                *DIST,
                "total_mass_before",
                "total_mass_after",
                "cell_coverage_before",
                "cell_coverage_after",
            ]
        ]
        .sort_values(["subset", "donor"])
        .reset_index(drop=True)
    )
    primary_a = a[a.subset.eq("all_T")]
    expected = data["dev_evidence"]["primary_multiscale"]["RFU"]
    require(
        np.allclose(
            primary_a[DIST].median(),
            [
                expected["d_clone_median"],
                expected["d_group_median"],
                expected["cancellation_median"],
            ],
        ),
        "Frozen GSE190905 medians differ.",
    )
    support = data["a_support"].query("policy == 'threshold'").reset_index(drop=True)
    require(
        support[support.subset.isin(["all_T", "CD4", "CD8"])][["before_cells", "after_cells"]]
        .ge(100)
        .all()
        .all(),
        "Changed support threshold.",
    )
    sensitivity = sensitivity_table(
        primary_a, data["a_distances"], data["a_depth"], data["a_dominance"]
    )
    bvis = data["b_visits"].copy()
    require(len(bvis) == 18 and bvis.source_available.sum() == 11, "External visit scope differs.")
    require(
        bvis.primary_TRB_cells.sum() == 13077 and bvis.qualified_cells.sum() == 9827,
        "External count totals differ.",
    )
    require(
        bvis.loc[~bvis.source_available, ["primary_TRB_cells", "qualified_cells", "coverage"]]
        .isna()
        .all()
        .all(),
        "Missing external visits were zero-imputed.",
    )
    require(
        np.allclose(bvis.qualified_cells / bvis.primary_TRB_cells, bvis.coverage, equal_nan=True),
        "External coverage denominator differs.",
    )
    b = primary(data["b_distances"])
    require(len(b) == 8, "External supported pairs differ.")
    pairs = external_pairs(data["b_pairs"], b)
    persistence = data["b_persistence"].copy()
    for num, den, ratio in [
        ("persistent_RFUs", "detected_union_RFUs", "persistent_over_union"),
        (
            "persistent_without_shared_receptor",
            "persistent_RFUs",
            "persistent_without_shared_receptor_fraction",
        ),
    ]:
        require(
            np.allclose(
                persistence[num] / persistence[den].replace(0, np.nan),
                persistence[ratio],
                equal_nan=True,
            ),
            "External persistence denominator differs.",
        )
    cross = data["cross"].copy()
    require(len(cross) == 7, "Cross-application frozen strata differ.")
    for row in cross.itertuples():
        if row.dataset == "RP1-14":
            val = next(
                r
                for r in data["rp_evidence"]["primary_by_compartment_grouping"]
                if r["compartment"] == row.compartment and r["grouping"] == "RFU"
            )
            measured = [val["clone_tv_median"], val["group_tv_median"], val["cancellation_median"]]
        else:
            part = (
                primary_a
                if row.dataset == "GSE190905"
                else b[b.compartment.eq(row.compartment) & b.interval.eq(row.interval)]
            )
            measured = part[DIST].median()
            require(
                part.donor.nunique() == row.participants, "Cross-application donor count differs."
            )
        require(
            np.allclose(
                measured,
                [row.receptor_TV_median, row.RFU_TV_median, row.cancellation_median],
                rtol=1e-12,
                atol=1e-12,
            ),
            "Cross-application endpoint differs.",
        )
    reg = data["reg_evidence"]
    counts = [
        reg["rfu_qtl_variants"],
        reg["eqtl"]["independent"]["variants"],
        reg["caqtl"]["independent"]["variants"],
        reg["both_independent"]["variants"],
    ]
    require(counts == [623, 17, 32, 6], "Frozen regulatory counts differ.")
    reg_counts = pd.DataFrame(
        {
            "category": [
                "published_RFU_QTL",
                "conditional_eQTL",
                "conditional_caQTL",
                "same_variant_both",
            ],
            "variants": counts,
            "denominator_variants": 623,
        }
    )
    matrix, targets = regulatory_matrix(data["reg_candidates"], data["reg_targets"])
    prediction, pred_summary = prediction_tables(data)
    w = data["w_coverage"]
    require(
        w.donor.nunique() == 24 and w.tissue.nunique() == 10, "Wells atlas coverage scope differs."
    )
    require(
        w.atlas_cells.sum() == 610429
        and w.primary_trb_cells.sum() == 303088
        and w.threshold_cells.sum() == 233913,
        "Wells counts differ.",
    )
    require(w.groupby("donor").primary_trb_cells.sum().eq(0).sum() == 3, "Zero-TRB donors lost.")
    wg = wells_grid(w)
    tissue = w.groupby("tissue", as_index=False).agg(
        atlas_cells=("atlas_cells", "sum"),
        primary_trb_cells=("primary_trb_cells", "sum"),
        threshold_cells=("threshold_cells", "sum"),
        atlas_donors=("donor", "nunique"),
    )
    support_w = data["w_support"].copy()
    require(
        support_w.public_support_rule.eq(support_w.cells.ge(50) & support_w.donors.ge(4)).all(),
        "Wells support rule changed.",
    )
    require(
        len(support_w) == 4928 and support_w.public_support_rule.sum() == 1471,
        "Wells RFU support differs.",
    )
    ctx = data["w_context"]
    require(
        ctx.cells.sum() == 233913 and ctx.cell_type.nunique() == 17,
        "Wells source annotation coverage differs.",
    )
    breadth = support_w.merge(
        ctx.groupby("rfu_label").cell_type.nunique().rename("cell_types"),
        on="rfu_label",
        validate="one_to_one",
    )
    require(
        np.array_equal(
            breadth.set_index("rfu_label").tissues.sort_index(),
            ctx.groupby("rfu_label").tissue.nunique().sort_index(),
        ),
        "RFU tissue support differs.",
    )
    histogram = pd.concat(
        [
            breadth.groupby(col)
            .size()
            .rename("rfus")
            .reset_index()
            .rename(columns={col: "categories_observed"})
            .assign(dimension=col)
            for col in ["tissues", "cell_types"]
        ],
        ignore_index=True,
    )
    celltypes = ctx.groupby("cell_type", as_index=False).agg(
        qualified_cells=("cells", "sum"),
        rfus=("rfu_label", "nunique"),
        donors=("donor", "nunique"),
        tissues=("tissue", "nunique"),
    )
    return {
        "A_visits": visits,
        "A_pairs": a,
        "A_support": support,
        "A_sensitivity": sensitivity,
        "A_subsamples": data["a_depth"][["donor", "replicate", "seed", "cells_per_visit", *DIST]],
        "B_visits": bvis,
        "B_pairs": pairs,
        "B_persistence": persistence,
        "C_summary": cross,
        "D_counts": reg_counts,
        "D_matrix": matrix,
        "D_target_records": targets,
        "E_all_weightings": prediction,
        "E_summary": pred_summary,
        "F_donor_tissue": wg,
        "F_tissues": tissue,
        "F_rfu_breadth": breadth,
        "F_breadth_histogram": histogram,
        "F_source_cell_types": celltypes,
    }


def write_metadata(out, tables, inputs, captions, command):
    schema = json.loads(SCHEMA.read_text())
    require(
        set().union(*(set(frame) for frame in tables.values())) <= schema.keys(),
        "Undocumented source column.",
    )
    json_save(out / "manifests/table_schema.json", schema)
    index = {}
    for key in "ABCDEF":
        names = [n for n in tables if n.startswith(key + "_")]
        sources = {n: inputs[n] for n in INPUTS_BY_FIGURE[key]}
        lines = [
            f"# {key}: {captions[key]['title']}",
            "",
            RULES[key],
            "",
            f"Generate: `{command}`",
            "",
            "## Exact completed inputs",
            "",
        ]
        lines += [f"- `{r['path']}` — SHA256 `{r['sha256']}`." for r in sources.values()]
        lines += [
            "",
            "## Numeric definitions",
            "",
            "`d_clone`/`receptor_TV_median` = Receptor TV (nucleotide CDR3/V/J); `d_group`/`rfu_tv`/`RFU_TV_median` = RFU TV. Cancellation is the within-donor receptor-minus-group difference. TV values are dimensionless [0,1]. Count units are explicit below; donors/technical replicates are never pooled. Blank numeric cells are undefined/unavailable, never imputed zero.",
            "",
            "Coverage fractions preserve their named denominator. Persistence is persistent_RFUs / detected_union_RFUs, and persistent_without_shared_receptor / persistent_RFUs. E uses clone-weighted log loss as primary and equal-donor means; positive delta is worse. D counts unique variants and deduplicated variant/layer/target/rank records, not joined RFU rows. F categories are unchanged source tissue/cell_type annotations. All other source fields retain their original names and units.",
            "",
        ]
        for name in names:
            frame = tables[name]
            lines += [
                f"## {name}.tsv",
                "",
                f"{len(frame)} rows. Columns (in stored order):",
                "",
                "| Column | Definition |",
                "| --- | --- |",
            ]
            lines += [f"| `{col}` | {schema[col]} |" for col in frame]
            lines += [""]
        lines += [
            "Figure-specific legend and take-home message:",
            "",
            captions[key]["manuscript_legend"],
            "",
            captions[key]["take_home"],
            "",
        ]
        (out / "tables" / f"{key}_README.md").write_text("\n".join(lines))
        index[key] = {
            **captions[key],
            "sources": sources,
            "source_tables": [f"tables/{n}.tsv" for n in names],
            "readme": f"tables/{key}_README.md",
            "rules": RULES[key],
            "command": command,
        }
    for key, components in [("poster_block", "ABD"), ("manuscript_assembly", "ABC")]:
        index[key] = {
            **captions[key],
            "components": list(components),
            "sources": {n: inputs[n] for k in components for n in INPUTS_BY_FIGURE[k]},
            "source_tables": sorted({n for k in components for n in index[k]["source_tables"]}),
            "readme": "README.md",
            "command": command,
        }
    lines = ["# Cross-application reusable captions", ""]
    for key, entry in captions.items():
        lines += [f"## {key} — {entry['title']}", ""]
        for field in ["poster_caption", "manuscript_legend", "take_home", "recommended_use"]:
            lines += [f"**{field.replace('_', ' ').capitalize()}:** {entry[field]}", ""]
    (out / "manifests/captions.md").write_text("\n".join(lines))
    json_save(out / "manifests/captions.json", captions)
    return index


def validate(out: Path, completion: dict, tables: dict) -> dict:
    from PIL import Image

    for name, checksum in completion["outputs"].items():
        path = safe_child(out, name)
        require(path.is_file() and sha256(path) == checksum, f"Output hash failed: {path}")
    for name, expected in tables.items():
        actual = pd.read_csv(out / "tables" / f"{name}.tsv", sep="\t")
        pd.testing.assert_frame_equal(
            actual,
            expected.reset_index(drop=True),
            check_dtype=False,
            check_exact=False,
            rtol=1e-12,
            atol=1e-12,
        )
    index = json.loads((out / "manifests/figure_index.json").read_text())
    for entry in index.values():
        for name in [entry["readme"], *entry["source_tables"]]:
            require(safe_child(out, name).is_file(), f"Missing figure source: {name}")
        for record in entry["sources"].values():
            path = Path(record["path"])
            require(
                path.is_file() and sha256(path) == record["sha256"],
                f"Changed input reference: {path}",
            )
    exports = json.loads((out / "manifests/exports.json").read_text())
    require(
        len(exports) == 16,
        "Expected twelve audience figures, two compact panels and two assemblies.",
    )
    require(
        {(r["figure"], r["mode"]) for r in exports if r["figure"] in "ABCDEF"}
        == {(key, mode) for key in "ABCDEF" for mode in ["manuscript", "poster"]},
        "Missing audience figure.",
    )
    for record in exports:
        require(
            not record["text_outside_canvas"] and not record["text_overlaps"],
            "Failed figure layout.",
        )
        require(
            record["minimum_font_pt"] >= (18 if record["mode"] == "poster" else 8),
            "Undersized export labels.",
        )
        for name in record["files"]:
            path = safe_child(out, name)
            require(name in completion["outputs"], f"Unhashed export: {name}")
            if path.suffix == ".pdf":
                require(path.read_bytes().startswith(b"%PDF-"), f"Invalid PDF: {path}")
            elif path.suffix == ".svg":
                root = ET.parse(path).getroot()
                require(
                    root.tag.endswith("svg")
                    and not any(el.tag.endswith("}image") for el in root.iter()),
                    f"SVG contains raster: {path}",
                )
            else:
                with Image.open(path) as img:
                    require(
                        img.size == tuple(round(v * 300) for v in record["size_inches"]),
                        f"PNG dimensions differ: {path}",
                    )
                    require(
                        all(abs(value - 300) < 0.1 for value in img.info.get("dpi", (0, 0))),
                        f"PNG dpi differs: {path}",
                    )
                    img.verify()
                if "biorender_blank_rectangle_fraction" in record:
                    left, bottom, width, height = record["biorender_blank_rectangle_fraction"]
                    with Image.open(path) as img:
                        box = (
                            round(left * img.width),
                            round((1 - bottom - height) * img.height),
                            round((left + width) * img.width),
                            round((1 - bottom) * img.height),
                        )
                        require(
                            img.crop(box).convert("RGB").getextrema() == ((255, 255),) * 3,
                            "Reserved BioRender region is not blank.",
                        )
    return {
        "status": "passed",
        "source_tables": len(tables),
        "export_sets": len(exports),
        "hashed_outputs": len(completion["outputs"]),
        "biological_analyses_rerun": False,
    }


def run(args) -> None:
    config, captions = json.loads(CONFIG.read_text()), json.loads(CAPTIONS.read_text())
    workspaces = {
        "methods": args.workspace.resolve(),
        "regulatory": args.regulatory_workspace.resolve(),
    }
    out = (args.out or workspaces["methods"] / "results" / config["figure_set"]).resolve()
    require(
        not out.is_relative_to(REPO), "Biological plotting tables/exports must remain outside Git."
    )
    for spec in config["upstream"].values():
        require(
            not out.is_relative_to(workspaces[spec["workspace"]] / spec["directory"]),
            "Cannot modify completed source directories.",
        )
    data, inputs, upstream = verify_sources(workspaces, config)
    tables = build_tables(data)
    paths = [
        Path(__file__),
        Path(__file__).with_name("_cross_application_drawing.py"),
        Path(__file__).with_name("_rp1_14_figure_drawing.py"),
        Path(__file__).with_name("rp1_14_reusable_figures.py"),
        CONFIG,
        CAPTIONS,
        SCHEMA,
    ]
    identity = {
        "inputs": inputs,
        "upstream": upstream,
        "code_sha256": {str(p.relative_to(REPO)): sha256(p) for p in paths},
        "versions": {
            n: importlib.metadata.version(n) for n in ["matplotlib", "pandas", "numpy", "Pillow"]
        },
        "python": sys.version,
        "executable": sys.executable,
    }
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    marker = out / "completion.json"
    if marker.exists():
        existing = json.loads(marker.read_text())
        require(
            existing["status"] == "complete" and existing["fingerprint"] == fingerprint,
            "Completed figure identity changed; use a new versioned directory.",
        )
        print(json.dumps(validate(out, existing, tables), sort_keys=True))
        print("Completed quantitative figures verified; no files rewritten.")
        return
    require(not args.verify_only, "No completed figure pack to verify.")
    require(
        not out.exists() or not any(out.iterdir()),
        "Output nonempty/incomplete; use a new directory.",
    )
    for name in [
        "figures/manuscript",
        "figures/poster",
        "figures/panels",
        "tables",
        "manifests",
        "preview",
    ]:
        (out / name).mkdir(parents=True, exist_ok=True)
    command = shlex.join(
        [
            sys.executable,
            "-m",
            "manuscript.scripts.cross_application_reusable_figures",
            "--workspace",
            str(workspaces["methods"]),
            "--regulatory-workspace",
            str(workspaces["regulatory"]),
            "--out",
            str(out),
        ]
    )
    for name, frame in tables.items():
        frame.to_csv(out / "tables" / f"{name}.tsv", sep="\t", index=False, float_format="%.17g")
    index = write_metadata(out, tables, inputs, captions, command)
    from ._cross_application_drawing import render_all

    exports = render_all(out, tables, config, captions)
    for key, entry in index.items():
        entry["exports"] = [r for r in exports if r["figure"] == key]
        json_save(out / "manifests" / f"{key}.json", entry)
    json_save(out / "manifests/figure_index.json", index)
    json_save(out / "manifests/exports.json", exports)
    json_save(
        out / "manifests/provenance.json",
        {
            **identity,
            "fingerprint": fingerprint,
            "command": command,
            "git_parent_at_generation": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
            ).strip(),
            "source_config": config,
            "new_biological_analyses": [],
            "conceptual_artwork": False,
        },
    )
    (out / "README.md").write_text(
        "# Cross-application reusable quantitative figures v1\n\n"
        "Data-only plots and summaries from completed scRFU outputs. No conceptual cartoons or causal arrows.\n\n"
        "A: GSE190905; B: GSE280982; C: cross-dataset transportability; D: Matos/RfuWAS; E: negative held-out prediction; F: Wells context.\n\n"
        "`figures/manuscript/` and `figures/poster/` contain each A–F figure in vector PDF/SVG and 300-dpi PNG. "
        "`figures/panels/` contains compact quantitative A/B components used in the manuscript assembly. "
        "`tables/` contains 19 TSVs plus per-figure source/definition READMEs. "
        "`manifests/` provides input/output hashes, figure selections, captions and take-home messages. "
        "`preview/` provides color/grayscale previews and contact sheets.\n\n"
        "Poster block: `figures/poster/cross_application_poster_block.pdf` (also SVG/PNG), 762 × 609.6 mm, minimum 18-pt text at native size. "
        "The lower-right normalized rectangle [0.63, 0.055, 0.35, 0.365] is deliberately blank for BioRender. "
        "The assembly manifest records this region; no placeholder artwork was drawn.\n\n"
        "Manuscript assembly: `figures/manuscript/cross_application_manuscript_assembly.pdf` (also SVG/PNG), "
        "190.5 × 266.7 mm, minimum 8-pt text. A/B/C are the main applications/synthesis; D/F are supporting context and E is supplementary/backup.\n\n"
        f"Regenerate from the repository root:\n\n```bash\n{command}\n```\n\n"
        "Add `--verify-only` to validate without rendering. Normal repeats also verify without rewriting. Changed inputs/code/software require a new output directory. "
        "All eight completed upstream manifests are verified, including every file in the protected RP1-14 figure pack. "
        "No endpoint, QTL search, assignment, prediction fit or random draw is rerun.\n\n"
        "Missing external visits remain missing, zero-TRB atlas samples remain in the denominator, technical draws are not donors, "
        "regulatory records are deduplicated across RFU joins, and read/cell effects are not pooled. "
        "TV = total variation; RFU = receptor functional unit; TRB = T-cell receptor beta; GEX = gene expression.\n"
    )
    initial = {
        "outputs": {
            str(p.relative_to(out)): sha256(p) for p in sorted(out.rglob("*")) if p.is_file()
        }
    }
    report = validate(out, initial, tables)
    # The final completion also hashes the validation report written below.
    report["hashed_outputs"] += 1
    report["upstream_outputs_verified"] = sum(v["verified_outputs"] for v in upstream.values())
    report["protected_rp1_14_pack_unchanged"] = True
    json_save(out / "manifests/validation.json", report)
    outputs = {str(p.relative_to(out)): sha256(p) for p in sorted(out.rglob("*")) if p.is_file()}
    json_save(marker, {"status": "complete", "fingerprint": fingerprint, "outputs": outputs})
    print(json.dumps(report, sort_keys=True))
    print(f"Reusable quantitative figures complete: {out}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--regulatory-workspace", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--verify-only", action="store_true")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
