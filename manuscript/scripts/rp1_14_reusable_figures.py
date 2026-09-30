"""Render reusable RP1-14 figures from hash-pinned results; never run analysis.

Only small frozen summary tables are parsed. No receptor sequences, assignments,
raw inputs, permutations, or subsampling routines are loaded or executed.
"""

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

REPO = Path(__file__).resolve().parents[2]
CONFIG = REPO / "manuscript/config/rp1_14_reusable_figures_v1.json"
CAPTIONS = REPO / "manuscript/figures/rp1_14_reusable_captions.json"
COMPARTMENTS = ("CD4", "CD8")
PAIR = ["donor", "compartment"]
ENDPOINTS = ["sample_before", "sample_after", "visit_before", "visit_after", "interval_years"]
METRICS = {
    "d_clone": "receptor_tv",
    "d_group": "rfu_tv",
    "aggregation_cancellation": "cancellation",
}
INPUTS = {
    "coverage": ("rp1_14_v1", "rp1_14_sample_coverage.tsv"),
    "multiscale": ("rp1_14_v1", "rp1_14_multiscale_pairs.tsv"),
    "controls": ("rp1_14_v1", "rp1_14_fixed_group_controls.tsv"),
    "persistence": ("rp1_14_v1", "rp1_14_persistence_donor_summary.tsv"),
    "sampling": ("rp1_14_v1", "rp1_14_empirical_read_sensitivity.tsv"),
    "evidence": ("rp1_14_v1", "evidence_counts.json"),
    "provenance": ("rp1_14_v1", "provenance.json"),
    "repair": ("rp1_14_repair_qc_v1", "evidence_counts.json"),
    "conflicts": ("rp1_14_repair_qc_v1", "repeated_receptor_conflicts.tsv"),
    "parity": ("rp1_14_repair_qc_v1", "parity_summary.tsv"),
}
FIGURE_INPUTS = {
    "A": ["coverage", "evidence", "provenance"],
    "B": ["multiscale", "evidence", "provenance"],
    "C": ["multiscale", "controls", "provenance"],
    "D": ["persistence", "provenance"],
    "E": ["multiscale", "sampling", "provenance"],
    "F": ["repair", "conflicts", "parity"],
}
FIGURE_RULES = {
    "A": "All 36 samples. Depths sum three visits within donor/compartment. Global unique counts come from frozen evidence_counts, never sums of sample richness.",
    "B": "primary_pair=True, policy=threshold, weighting=reads. Same qualified universe for RFU, TRBV and TRBV_TRBJ. Equal-donor medians; n=6 per compartment.",
    "C": "Same primary RFU rows as B. Each donor's 30 already-frozen matched maps summarized with linear 0.025/0.5/0.975 quantiles; no new permutations or p-values.",
    "D": "Frozen donor summary for primary pairs, threshold reads, detection >=5 reads in BOTH endpoints. Numerator=no shared observed nucleotide-CDR3/V/J identity; denominator=persistent RFUs. Undefined fractions remain undefined.",
    "E": "Primary threshold RFU rows with reads or unique_clone weights. 50 existing observed-read subsamples per pair, 500 reads/visit. First summarize within donor, then across six donors; do not pool replicates as donors.",
    "F": "Completed repair evidence/conflict/parity summaries only. Near threshold is [0.59,0.61). Factor strata partition the same 1600 assays; never sum across factors.",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_save(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def safe_child(root: Path, name: str) -> Path:
    path = (root / name).resolve()
    require(path.is_relative_to(root.resolve()), f"Path escapes manifest root: {name}")
    return path


def verify_upstream(results: Path, config: dict) -> dict:
    """Verify pinned manifests AND every recorded output, including unused ones."""
    verified = {}
    for stage, expected in config["upstream"].items():
        root = results / stage
        manifest = root / "completion.json"
        require(manifest.is_file(), f"Missing upstream completion: {manifest}")
        require(sha256(manifest) == expected, f"Frozen completion hash changed: {manifest}")
        record = json.loads(manifest.read_text())
        require(record.get("status") == "complete", f"Upstream not complete: {stage}")
        for name, checksum in record["outputs"].items():
            path = safe_child(root, name)
            require(path.is_file() and sha256(path) == checksum, f"Upstream hash failed: {path}")
        verified[stage] = {
            "path": str(manifest),
            "sha256": expected,
            "verified_outputs": len(record["outputs"]),
        }
    return verified


def load_inputs(results: Path) -> tuple[dict, dict]:
    data, inputs = {}, {}
    for key, (stage, name) in INPUTS.items():
        path = results / stage / name
        data[key] = (
            json.loads(path.read_text()) if path.suffix == ".json" else pd.read_csv(path, sep="\t")
        )
        inputs[key] = {"path": str(path), "sha256": sha256(path), "upstream": stage}
    return data, inputs


def primary_rows(table: pd.DataFrame, *, weighting: str = "reads") -> pd.DataFrame:
    require(pd.api.types.is_bool_dtype(table.primary_pair), "primary_pair must be boolean.")
    out = table[table.primary_pair & table.policy.eq("threshold") & table.weighting.eq(weighting)]
    require(out.status.eq("valid").all(), "Invalid frozen primary pair.")
    require(
        out.visit_before.eq(1).all() and out.visit_after.eq(3).all(), "Changed primary endpoints."
    )
    require(not out.duplicated(PAIR + ["grouping"]).any(), "Duplicate primary measurement.")
    return out.copy()


def verify_pairs(frame: pd.DataFrame, donors: list[str], replicates: int = 1) -> None:
    expected = pd.MultiIndex.from_product([donors, COMPARTMENTS], names=PAIR)
    counts = frame.groupby(PAIR).size().reindex(expected)
    require(
        len(frame) == len(expected) * replicates and counts.eq(replicates).all(),
        f"Incomplete or duplicated donor/compartment coverage (expected {replicates} rows/pair).",
    )


def summarize_donors(frame: pd.DataFrame, metrics: list[str], by: list[str]) -> pd.DataFrame:
    """Equal-donor descriptive summaries; never pool resampling replicates."""
    rows = []
    for key, group in frame.groupby(by, sort=True):
        key = key if isinstance(key, tuple) else (key,)
        require(not group.donor.duplicated().any(), "Summarize within donor before across donors.")
        for metric in metrics:
            rows.append(
                {
                    **dict(zip(by, key, strict=True)),
                    "metric": metric,
                    "n_donors": group.donor.nunique(),
                    "median": group[metric].median(),
                    "minimum": group[metric].min(),
                    "maximum": group[metric].max(),
                }
            )
    return pd.DataFrame(rows)


def replicate_summary(frame: pd.DataFrame, donors: list[str], n: int) -> pd.DataFrame:
    verify_pairs(frame, donors, n)
    require(not frame.duplicated(PAIR + ["replicate"]).any(), "Duplicate frozen replicate.")
    require(
        frame.primary_pair.all() and frame.status.eq("valid").all(), "Changed sensitivity scope."
    )
    require(
        all(set(g.replicate) == set(range(n)) for _, g in frame.groupby(PAIR)),
        "Frozen replicate indices changed.",
    )
    rows = []
    for (donor, comp), group in frame.groupby(PAIR, sort=True):
        row = {"donor": donor, "compartment": comp, "n_replicates": n}
        for old, new in METRICS.items():
            lo, med, hi = group[old].quantile([0.025, 0.5, 0.975], interpolation="linear")
            row.update({new: med, f"{new}_lo": lo, f"{new}_hi": hi})
        rows.append(row)
    return pd.DataFrame(rows)


def build_tables(data: dict, config: dict) -> dict[str, pd.DataFrame]:
    """Figure-specific projections and descriptive summaries of frozen rows."""
    require(
        config["selection"]
        == {
            "primary_pair": True,
            "policy": "threshold",
            "weighting": "reads",
            "threshold": 0.6,
            "persistence_detection_reads": 5,
            "control_replicates": 30,
            "subsample_replicates": 50,
            "subsample_reads_per_visit": 500,
            "quantiles": [0.025, 0.5, 0.975],
            "quantile_interpolation": "linear",
        },
        "Figure selection is frozen; this runner cannot change endpoint definitions.",
    )
    coverage, evidence, multiscale = data["coverage"], data["evidence"], data["multiscale"]
    for key, expected in config["expected"].items():
        if key != "visits_per_donor_compartment":
            require(evidence[key] == expected, f"Frozen evidence count changed: {key}")
    frozen = data["provenance"]["config"]
    require(
        frozen["threshold"] == 0.6
        and frozen["rfu_detection_reads"] == [1, 5]
        and frozen["random_group_replicates"] == 30
        and frozen["depth_replicates"] == 50
        and frozen["depth_cap_reads"] == 500,
        "Frozen definitions differ from the figure contract.",
    )
    donors = sorted(coverage.donor.unique())
    require(len(donors) == 6 and len(coverage) == 36, "RP1-14 design changed.")
    require(not coverage.duplicated([*PAIR, "visit"]).any(), "Duplicate visit.")
    verify_pairs(coverage, donors, 3)
    require(
        all(set(g.visit) == {1, 2, 3} for _, g in coverage.groupby(PAIR)),
        "Incomplete chronological visits.",
    )
    timeline = coverage.pivot(
        index=["donor", "visit"], columns="compartment", values="elapsed_years"
    )
    require(np.allclose(timeline.CD4, timeline.CD8), "Compartment visit dates disagree.")
    samples = (
        coverage[
            [
                "sample_id",
                *PAIR,
                "visit",
                "elapsed_years",
                "all_source_reads",
                "productive_read_depth",
                "fixed_primary_reads",
                "fixed_primary_receptors",
                "fixed_primary_rfus",
                "fixed_primary_fraction_all_reads",
                "fixed_primary_fraction_productive_reads",
            ]
        ]
        .sort_values([*PAIR, "visit"])
        .reset_index(drop=True)
    )
    for denom, col in [
        ("all_source_reads", "fixed_primary_fraction_all_reads"),
        ("productive_read_depth", "fixed_primary_fraction_productive_reads"),
    ]:
        require(
            np.allclose(samples.fixed_primary_reads / samples[denom], samples[col]),
            "Coverage denominator differs.",
        )
    require(
        samples.all_source_reads.sum() == evidence["source_reads"], "Source read total differs."
    )
    require(
        samples.fixed_primary_reads.sum() == evidence["primary_mapped_reads"],
        "Qualified read total differs.",
    )
    depth = samples.groupby(PAIR, as_index=False)[["all_source_reads", "fixed_primary_reads"]].sum()
    depth["qualified_fraction_source"] = depth.fixed_primary_reads / depth.all_source_reads
    totals = pd.DataFrame(
        [
            {"metric": key, "value": evidence[key]}
            for key in config["expected"]
            if key != "visits_per_donor_compartment"
        ]
        + [
            {"metric": f"follow_up_years_{s}", "value": evidence["interval_years_min_max"][i]}
            for i, s in enumerate(["min", "max"])
        ]
    )
    primary = primary_rows(multiscale)
    require(set(primary.grouping) == {"RFU", "TRBV", "TRBV_TRBJ"}, "Grouping definitions differ.")
    rfu = primary[primary.grouping.eq("RFU")].sort_values(PAIR)
    verify_pairs(rfu, donors)
    pairs = rfu[[*PAIR, *ENDPOINTS, *METRICS]].rename(columns=METRICS).reset_index(drop=True)
    for grouping, label in [("TRBV", "trbv_tv"), ("TRBV_TRBJ", "vj_tv")]:
        other = primary[primary.grouping.eq(grouping)]
        verify_pairs(other, donors)
        require(
            np.allclose(
                other.set_index(PAIR).d_clone.sort_index(), rfu.set_index(PAIR).d_clone.sort_index()
            ),
            "Representations do not share receptor TV.",
        )
        pairs = pairs.merge(
            other[[*PAIR, "d_group"]].rename(columns={"d_group": label}),
            on=PAIR,
            validate="one_to_one",
        )
    require(
        np.allclose(pairs.receptor_tv - pairs.rfu_tv, pairs.cancellation),
        "Cancellation identity failed.",
    )
    b_summary = summarize_donors(
        pairs, ["receptor_tv", "rfu_tv", "trbv_tv", "vj_tv", "cancellation"], ["compartment"]
    )
    for expected in evidence["primary_by_compartment_grouping"]:
        match = primary[
            primary.compartment.eq(expected["compartment"])
            & primary.grouping.eq(expected["grouping"])
        ]
        for col, field in [
            ("d_clone", "clone_tv_median"),
            ("d_group", "group_tv_median"),
            ("aggregation_cancellation", "cancellation_median"),
        ]:
            require(
                np.isclose(match[col].median(), expected[field], rtol=1e-12, atol=1e-12),
                f"Frozen median mismatch: {field}",
            )
    controls = data["controls"]
    control_summary = replicate_summary(controls, donors, 30)
    c_donors = (
        pairs[[*PAIR, "cancellation"]]
        .rename(columns={"cancellation": "observed_cancellation"})
        .merge(
            control_summary[
                [*PAIR, "n_replicates", "cancellation", "cancellation_lo", "cancellation_hi"]
            ].rename(
                columns={
                    "cancellation": "control_median",
                    "cancellation_lo": "control_q025",
                    "cancellation_hi": "control_q975",
                }
            ),
            on=PAIR,
            validate="one_to_one",
        )
    )
    control_base = controls.groupby(PAIR).d_clone.first().sort_index()
    require(
        np.allclose(control_base, rfu.set_index(PAIR).d_clone.sort_index()),
        "Control receptor TV differs.",
    )
    persistence = data["persistence"].sort_values(PAIR).reset_index(drop=True).copy()
    verify_pairs(persistence, donors)
    for numerator, denom, frac in [
        ("persistent_rfus", "union_detected_rfus", "persistent_fraction_union"),
        (
            "persistent_without_shared_clone",
            "persistent_rfus",
            "fraction_persistent_without_shared_clone",
        ),
    ]:
        calculated = persistence[numerator].div(persistence[denom].replace(0, np.nan))
        require(
            np.allclose(calculated, persistence[frac], equal_nan=True, rtol=1e-12),
            "Persistence denominator differs.",
        )
    persistence["detection_reads_per_endpoint"] = 5
    d_summary = summarize_donors(
        persistence,
        ["persistent_fraction_union", "fraction_persistent_without_shared_clone"],
        ["compartment"],
    )
    schematic = pd.DataFrame(
        [
            {
                "hypothetical_rfu": "same_group",
                "visit": visit,
                "arbitrary_receptor_label": label,
                "is_observed_data": False,
            }
            for visit, labels in [("early", "ab"), ("late", "cd")]
            for label in labels
        ]
    )
    sampling = data["sampling"]
    require(sampling.reads_per_visit.eq(500).all(), "Subsample depth differs.")
    require(
        sampling.scheme.eq(
            "observed_reads_without_replacement_conditional_on_reuse_universe"
        ).all(),
        "Subsampling scheme differs.",
    )
    sensitivity = []
    for weighting, condition in [("reads", "read_weighted"), ("unique_clone", "unique_receptor")]:
        selected = primary_rows(multiscale, weighting=weighting)
        selected = selected[selected.grouping.eq("RFU")]
        verify_pairs(selected, donors)
        part = selected[[*PAIR, *METRICS]].rename(columns=METRICS).copy()
        part["condition"], part["n_replicates"] = condition, 1
        for metric in METRICS.values():
            part[f"{metric}_lo"] = part[metric]
            part[f"{metric}_hi"] = part[metric]
        sensitivity.append(part)
    sensitivity.append(replicate_summary(sampling, donors, 50).assign(condition="500_reads"))
    e_donors = (
        pd.concat(sensitivity, ignore_index=True)
        .sort_values([*PAIR, "condition"])
        .reset_index(drop=True)
    )
    e_summary = summarize_donors(e_donors, list(METRICS.values()), ["compartment", "condition"])
    repair, conflicts, parity = data["repair"], data["conflicts"].set_index("value"), data["parity"]
    all_parity = parity[parity.factor.eq("all") & parity.level.eq("all")].iloc[0]
    near = parity[
        parity.factor.eq("score_bin") & parity.level.isin(["0.59_to_0.60", "0.60_to_0.61"])
    ].n_unique_receptors.sum()
    require(repair["RFU_and_score_vectors_identical"], "Repair changed saved vectors.")
    require(
        all_parity.n_unique_receptors == repair["n_assayed_unique_receptors"] == 1600,
        "Parity assay size differs.",
    )
    require(near == 485, "Near-threshold assay count differs.")
    checks = [
        (
            "vector_entries_unchanged",
            repair["source_rows"],
            repair["source_rows"],
            "RFU and score vectors identical",
        ),
        (
            "amino_acid_labels_repaired",
            repair["AA_labels_repaired"],
            repair["source_rows"],
            "Row-label repair only",
        ),
        (
            "amino_acid_labels_unchanged",
            repair["AA_labels_unchanged"],
            repair["source_rows"],
            "Original labels retained",
        ),
        (
            "rfu_conflicts_before",
            conflicts.loc["RFU", "before_conflicting_AA_groups"],
            np.nan,
            "Distinct conflicting amino-acid groups",
        ),
        (
            "rfu_conflicts_after",
            conflicts.loc["RFU", "after_conflicting_AA_groups"],
            np.nan,
            "Distinct conflicting amino-acid groups",
        ),
        (
            "score_conflicts_before",
            conflicts.loc["score", "before_conflicting_AA_groups"],
            np.nan,
            "Distinct conflicting amino-acid groups",
        ),
        (
            "score_conflicts_after",
            conflicts.loc["score", "after_conflicting_AA_groups"],
            np.nan,
            "Distinct conflicting amino-acid groups",
        ),
        ("rfu_label_matches", 1600 - repair["label_mismatches"], 1600, "Frozen stratified parity"),
        (
            "threshold_matches",
            1600 - repair["threshold_mismatches"],
            1600,
            "Frozen score >=0.6 decisions",
        ),
        ("near_threshold_assays", near, 1600, "Historical score in [0.59,0.61)"),
        (
            "max_abs_score_difference",
            repair["max_abs_score_difference"],
            np.nan,
            "Absolute floating-point discrepancy",
        ),
        ("score_tolerance", repair["score_comparison_tolerance"], np.nan, "Frozen assay tolerance"),
    ]
    replicate_cols = [*PAIR, *ENDPOINTS, "replicate", "seed", *METRICS]
    return {
        "A_samples": samples,
        "A_donor_depth": depth,
        "A_totals": totals,
        "B_pairs": pairs,
        "B_summary": b_summary,
        "C_donors": c_donors,
        "C_maps": controls[
            replicate_cols + ["exchangeable_receptor_fraction", "changed_label_fraction"]
        ].rename(columns=METRICS),
        "D_donors": persistence,
        "D_summary": d_summary,
        "D_schematic": schematic,
        "E_donors": e_donors,
        "E_summary": e_summary,
        "E_subsamples": sampling[replicate_cols + ["reads_per_visit", "scheme"]].rename(
            columns=METRICS
        ),
        "F_checks": pd.DataFrame(checks, columns=["check", "value", "denominator", "definition"]),
        "F_parity_strata": parity,
    }


COLUMN_NOTES = {
    "donor": "Frozen RP01–RP06 label; the biological unit.",
    "compartment": "CD4 or CD8; the same six donors contribute both.",
    "sample_id": "Frozen donor/compartment/chronological-visit identifier.",
    "visit": "Chronological visit (1, 2, 3); schematic rows instead use early/late.",
    "elapsed_years": "Years from the donor's first collection.",
    "all_source_reads": "All source sequencing reads; summed over three visits in A_donor_depth only.",
    "productive_read_depth": "Productive sequencing-read denominator for this sample.",
    "fixed_primary_reads": "Qualified sequencing reads in the frozen reused dictionary, score >=0.6.",
    "fixed_primary_receptors": "Distinct qualified nucleotide-CDR3/V/J identities in this sample.",
    "fixed_primary_rfus": "Distinct qualified receptor functional units in this sample.",
    "fixed_primary_fraction_all_reads": "fixed_primary_reads / all_source_reads.",
    "fixed_primary_fraction_productive_reads": "fixed_primary_reads / productive_read_depth.",
    "qualified_fraction_source": "Ratio of donor/compartment summed qualified to summed source reads.",
    "sample_before": "Frozen earliest endpoint sample identifier.",
    "sample_after": "Frozen latest endpoint sample identifier.",
    "visit_before": "Chronological earliest endpoint, 1.",
    "visit_after": "Chronological latest endpoint, 3.",
    "interval_years": "Frozen elapsed years between endpoints.",
    "receptor_tv": "Receptor TV (nucleotide CDR3/V/J), dimensionless [0,1].",
    "rfu_tv": "RFU TV (or matched-group TV in C_maps), dimensionless [0,1].",
    "trbv_tv": "TRBV TV, using recorded variable-call categories including family/unresolved calls.",
    "vj_tv": "V/J TV, using recorded joint variable/joining-call categories.",
    "cancellation": "Aggregation cancellation: receptor_tv minus group TV; never a difference of compartment medians.",
    "observed_cancellation": "Frozen observed receptor_tv minus rfu_tv for one donor/compartment.",
    "control_median": "Within-donor median cancellation across 30 fixed matched maps.",
    "control_q025": "2.5th percentile of 30 matched-map cancellations, linear interpolation.",
    "control_q975": "97.5th percentile of 30 matched-map cancellations, linear interpolation.",
    "replicate": "Frozen zero-based technical replicate, not a donor.",
    "seed": "Already-used upstream random seed, retained for traceability; no new draws.",
    "exchangeable_receptor_fraction": "Fraction in frozen feature strata with more than one RFU label.",
    "changed_label_fraction": "Fraction of frozen labels changed by this matched map.",
    "n_replicates": "30 maps, 50 read subsamples, or 1 observed estimate, as applicable.",
    "persistent_rfus": "RFUs detected at >=5 reads in each endpoint.",
    "union_detected_rfus": "RFUs detected at >=5 reads at either endpoint.",
    "persistent_fraction_union": "persistent_rfus / union_detected_rfus.",
    "persistent_without_shared_clone": "Persistent RFUs with zero shared observed nucleotide-CDR3/V/J identities.",
    "fraction_persistent_without_shared_clone": "persistent_without_shared_clone / persistent_rfus; undefined if denominator is zero.",
    "median_dominant_fraction_before": "Frozen median within-RFU dominant-receptor read fraction at earliest endpoint.",
    "median_dominant_fraction_after": "Frozen median within-RFU dominant-receptor read fraction at latest endpoint.",
    "detection_reads_per_endpoint": "Unchanged RFU detection threshold of 5 qualified reads.",
    "hypothetical_rfu": "Same arbitrary group in the explanatory schematic only.",
    "arbitrary_receptor_label": "Invented a/b/c/d labels; not sequences or measured frequencies.",
    "is_observed_data": "False for every schematic row; no schematic quantitative inference.",
    "condition": "read_weighted, unique_receptor (one per observed identity/visit), or 500_reads.",
    "reads_per_visit": "Common frozen sensitivity depth, 500 observed reads per endpoint.",
    "scheme": "Frozen without-replacement observed-read sampling, conditional on reuse universe.",
    "metric": "Name of the metric summarized; see metric definitions here or A_totals frozen evidence keys.",
    "value": "Frozen numeric check/count; A_totals counts are global, not pooled sample richness.",
    "n_donors": "Six equally weighted biological donors; not number of technical replicates.",
    "median": "Equal-donor median; sensitivity first takes the median within each donor.",
    "minimum": "Smallest donor estimate; not a confidence bound.",
    "maximum": "Largest donor estimate; not a confidence bound.",
    "check": "Named frozen repair/parity check.",
    "denominator": "Applicable count denominator; blank when a ratio is not defined.",
    "definition": "Meaning and units of this repair check.",
    "factor": "Stratification factor partitioning the same parity assay.",
    "level": "Level within a single parity factor; never sum different factors.",
    "n_unique_receptors": "Unique amino-acid receptors in this assay stratum.",
    "label_mismatches": "RFU-label disagreements in the frozen assay.",
    "threshold_mismatches": "Disagreements in score >=0.6 decisions.",
    "score_mismatches_at_1e_12": "Absolute score discrepancies above 1e-12.",
    "max_abs_score_difference": "Maximum absolute score difference in this stratum.",
}


def column_note(column: str) -> str:
    if column.endswith(("_lo", "_hi")):
        return "Within-donor 2.5th/97.5th technical percentile (linear); equals estimate for unsampled conditions."
    require(column in COLUMN_NOTES, f"Undocumented source column: {column}")
    return COLUMN_NOTES[column]


def write_documentation(
    out: Path, tables: dict, inputs: dict, captions: dict, command: str
) -> dict:
    records = {}
    for key in "ABCDEF":
        sources = {name: inputs[name] for name in FIGURE_INPUTS[key]}
        names = [name for name in tables if name.startswith(key + "_")]
        body = [
            f"# Figure {key}: source tables",
            "",
            FIGURE_RULES[key],
            "",
            f"Generate: `{command}`",
            "",
            "## Exact frozen inputs",
            "",
        ]
        for name, record in sources.items():
            body.append(f"- `{record['path']}` — SHA256 `{record['sha256']}` ({name}).")
        for name in names:
            frame = tables[name]
            body.extend(
                [
                    "",
                    f"## {name}.tsv ({len(frame)} rows)",
                    "",
                    "| Column | Definition |",
                    "| --- | --- |",
                ]
            )
            body.extend(f"| `{col}` | {column_note(col)} |" for col in frame)
        body.extend(
            [
                "",
                "Empty numeric values are undefined, never zero-imputed. Source tables are external artifacts, not Git inputs.",
                "",
            ]
        )
        (out / "tables" / f"{key}_README.md").write_text("\n".join(body))
        records[key] = {
            **captions[key],
            "sources": sources,
            "selection_and_summary": FIGURE_RULES[key],
            "source_tables": [f"tables/{name}.tsv" for name in names],
            "table_readme": f"tables/{key}_README.md",
            "command": command,
        }
    caption_lines = [
        "# RP1-14 reusable figure captions",
        "",
        "Definitions and wording accompany the frozen figure exports.",
        "",
    ]
    for key, entry in captions.items():
        caption_lines += [f"## {key} — {entry['title']}", ""]
        for field in ["poster_caption", "manuscript_legend", "take_home", "recommended_use"]:
            caption_lines += [f"**{field.replace('_', ' ').capitalize()}:** {entry[field]}", ""]
    (out / "manifests/rp1_14_reusable_figure_captions.md").write_text("\n".join(caption_lines))
    json_save(out / "manifests/rp1_14_reusable_figure_captions.json", captions)
    return records


def validate_outputs(out: Path, manifest: dict, tables: dict | None = None) -> dict:
    """Check hashes, table values, figure formats and every manifest reference."""
    from PIL import Image

    checked = 0
    for name, checksum in manifest["outputs"].items():
        path = safe_child(out, name)
        require(path.is_file() and sha256(path) == checksum, f"Output hash failed: {path}")
        checked += 1
    figure_index = json.loads((out / "manifests/figure_index.json").read_text())
    for entry in figure_index.values():
        for name in entry["source_tables"] + [entry["table_readme"]]:
            require(safe_child(out, name).is_file(), f"Missing figure source: {name}")
        for record in entry["sources"].values():
            path = Path(record["path"])
            require(
                path.is_file() and sha256(path) == record["sha256"],
                f"Input reference failed: {path}",
            )
    exports = json.loads((out / "manifests/exports.json").read_text())
    for record in exports:
        for name in record["files"]:
            path = safe_child(out, name)
            require(name in manifest["outputs"], f"Export missing from output hashes: {name}")
            if path.suffix == ".pdf":
                require(path.read_bytes().startswith(b"%PDF-"), f"Bad PDF: {path}")
            elif path.suffix == ".svg":
                root = ET.parse(path).getroot()
                require(root.tag.endswith("svg"), f"Bad SVG: {path}")
                require(
                    not any(e.tag.endswith("}image") for e in root.iter()),
                    f"Raster content in SVG: {path}",
                )
            elif path.suffix == ".png":
                with Image.open(path) as img:
                    expected = tuple(round(v * record["dpi"]) for v in record["size_inches"])
                    require(
                        img.size == expected,
                        f"PNG dimensions differ: {path}: {img.size} != {expected}",
                    )
                    img.verify()
        require(
            record["minimum_font_pt"] >= record["required_minimum_font_pt"],
            "Font below declared minimum.",
        )
        require(not record["text_outside_canvas"], "Text outside figure canvas.")
        require(not record["text_overlaps"], "Overlapping figure text.")
    if tables is not None:
        for name, expected in tables.items():
            actual = pd.read_csv(out / "tables" / f"{name}.tsv", sep="\t")
            pd.testing.assert_frame_equal(
                actual, expected, check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-12
            )
    return {
        "status": "passed",
        "hashed_outputs_checked": checked,
        "export_sets": len(exports),
        "source_tables": len(tables) if tables else 15,
        "upstream_endpoints_recomputed": False,
    }


def run(args: argparse.Namespace) -> None:
    config = json.loads(CONFIG.read_text())
    out = (args.out or args.workspace / "results" / config["figure_set"]).resolve()
    results = args.workspace.resolve() / "results"
    require(not out.is_relative_to(REPO), "Figure/source-table outputs must remain outside Git.")
    require(
        not any(
            out == results / stage or out.is_relative_to(results / stage)
            for stage in config["upstream"]
        ),
        "Cannot write into frozen upstream results.",
    )
    upstream = verify_upstream(results, config)
    data, inputs = load_inputs(results)
    tables = build_tables(data, config)
    captions = json.loads(CAPTIONS.read_text())
    code_paths = [
        Path(__file__),
        Path(__file__).with_name("_rp1_14_figure_drawing.py"),
        CONFIG,
        CAPTIONS,
    ]
    identity = {
        "code_sha256": {str(p.relative_to(REPO)): sha256(p) for p in code_paths},
        "upstream": upstream,
        "inputs": inputs,
        "versions": {
            name: importlib.metadata.version(name)
            for name in ["numpy", "pandas", "matplotlib", "Pillow"]
        },
        "python": sys.version,
        "executable": sys.executable,
    }
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    completion = out / "completion.json"
    if completion.exists():
        existing = json.loads(completion.read_text())
        require(
            existing.get("status") == "complete" and existing["fingerprint"] == fingerprint,
            "Completed figure identity changed; use a new versioned output directory.",
        )
        print(json.dumps(validate_outputs(out, existing, tables), sort_keys=True))
        print(
            "Completed reusable figures verified; no files rewritten and no biological analysis run."
        )
        return
    require(not args.verify_only, "No completed reusable figures to verify.")
    require(
        not out.exists() or not any(out.iterdir()),
        "Output is nonempty/incomplete; use a new versioned directory.",
    )
    for directory in [
        "figures/manuscript",
        "figures/poster",
        "figures/panels",
        "tables",
        "manifests",
        "preview",
    ]:
        (out / directory).mkdir(parents=True, exist_ok=True)
    command = shlex.join(
        [
            sys.executable,
            "-m",
            "manuscript.scripts.rp1_14_reusable_figures",
            "--workspace",
            str(args.workspace.resolve()),
            "--out",
            str(out),
        ]
    )
    for name, frame in tables.items():
        frame.to_csv(out / "tables" / f"{name}.tsv", sep="\t", index=False, float_format="%.17g")
    figure_index = write_documentation(out, tables, inputs, captions, command)
    from ._rp1_14_figure_drawing import render_all

    exports = render_all(out, tables, config, captions)
    for key, components in [("poster_block", "ACDE"), ("manuscript_main", "ABCD")]:
        figure_index[key] = {
            **captions[key],
            "component_figures": list(components),
            "command": command,
            "sources": {name: inputs[name] for k in components for name in FIGURE_INPUTS[k]},
            "source_tables": sorted(
                {name for k in components for name in figure_index[k]["source_tables"]}
            ),
            "table_readme": "tables/README.md",
        }
    for key, entry in figure_index.items():
        entry["exports"] = [record for record in exports if record["figure"] == key]
        json_save(out / "manifests" / f"{key}.json", entry)
    json_save(out / "manifests/figure_index.json", figure_index)
    json_save(out / "manifests/exports.json", exports)
    json_save(
        out / "manifests/provenance.json",
        {
            **identity,
            "fingerprint": fingerprint,
            "command": command,
            "working_directory": str(REPO),
            "git_parent_at_generation": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
            ).strip(),
            "frozen_config": config,
            "new_biological_analyses": [],
            "rendered_from": "Figure-specific projections/descriptive summaries of frozen output tables only.",
        },
    )
    (out / "tables/README.md").write_text(
        "# Source tables\n\nEach A–F README defines every column, input hash, denominator and command.\n"
        "All tables are TSV, with 17 significant digits for numeric round trips. Blank values remain undefined.\n"
        "No receptor sequences or new biological measurements are included. D_schematic is explicitly hypothetical.\n"
        "C_maps and E_subsamples retain frozen replicate indices and seeds; summaries never count those as donors.\n"
    )
    (out / "README.md").write_text(
        "# Reusable RP1-14 figures v1\n\n"
        "Frozen source checkpoint: `9078dc3e9e6e7adb95e7cc675bc5cd5e014f242b`.\n\n"
        "A–F each have manuscript and poster PDF/SVG/300-dpi PNG exports. "
        "The PDF and SVG are vector drawings; SVG text is outlined for portable typography.\n\n"
        "- `figures/manuscript/`: A–F and the A–D main-figure draft (190.5 × 264.2 mm).\n"
        "- `figures/poster/`: A–F and the A/C/D/E poster block (762 × 508 mm).\n"
        "- `figures/panels/`: separately reusable programmatic schematic and compartment key.\n"
        "- `tables/`: source TSVs and one column/denominator README per figure.\n"
        "- `manifests/`: exact input paths/hashes, selections, captions, output sizes and font checks.\n"
        "- `preview/`: lightweight color/grayscale previews and contact sheets.\n\n"
        "Poster text is at least 18 pt at native export size; keep the block at least 762 mm wide "
        "for that minimum on A0. Scaling smaller also scales text; use individual panels to rearrange. "
        "Standalone poster panels are 355.6 × 223.5 mm. Manuscript text is at least 8 pt at native size.\n\n"
        "Use C and D for the central poster message, A for context and E as an inset. "
        "A–D form the main-paper draft; E and F remain supplementary, with F available for Q&A. "
        "The original six-panel `rp1_14_v1/rp1_14_longitudinal.pdf` and all upstream files are preserved unchanged. "
        "Its additional within/between-donor and paired-compartment views remain available there.\n\n"
        f"Regenerate from the repository root:\n\n```bash\n{command}\n```\n\n"
        "Repeat with `--verify-only` to check all sources and exports without writing files. "
        "Changed code, inputs or software require a new output version; completed exports are never overwritten.\n\n"
        "No completed assignment, biological endpoint, random map, or read draw was rerun. "
        "Persistence requires >=5 qualified reads per endpoint. Matched-map and read-subsampling "
        "ranges are technical variation, not donor confidence intervals. RFU persistence does not "
        "establish preserved antigen function.\n"
    )
    preliminary = {
        "outputs": {
            str(p.relative_to(out)): sha256(p) for p in sorted(out.rglob("*")) if p.is_file()
        }
    }
    validation = validate_outputs(out, preliminary, tables)
    validation["upstream_outputs_verified"] = sum(x["verified_outputs"] for x in upstream.values())
    validation["frozen_primary_summary_medians_checked"] = 18
    validation["source_round_trip_checked"] = True
    json_save(out / "manifests/validation.json", validation)
    outputs = {str(p.relative_to(out)): sha256(p) for p in sorted(out.rglob("*")) if p.is_file()}
    json_save(completion, {"status": "complete", "fingerprint": fingerprint, "outputs": outputs})
    print(json.dumps(validation, sort_keys=True))
    print(f"Reusable RP1-14 figures complete: {out}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--verify-only", action="store_true")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
