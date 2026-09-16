"""Reuse pinned public assignments for the radiotherapy methods manuscript.

Run as ``python -m manuscript.scripts.radiation_methods --help`` from the
repository. This example never assigns RFUs or accepts an unpublished cohort.
Inputs, cell identifiers, donor crosswalks, results and figures stay outside Git.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.metadata
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import scrfu
from scrfu.io import file_sha256
from scrfu.tl import (
    longitudinal_similarity,
    multiscale_repertoire_change,
    permute_fixed_groups,
    validate_completed_rfu_run,
    validate_longitudinal_design,
)

from ._figure_common import clean_axis, panel_label, save_figure, set_style

REPO = Path(__file__).resolve().parents[2]
ALLOWED = {
    "RP1-14_published_authorized",
    "Wells_public",
    "GSE190905",
    "GSE280982",
    "Matos_RfuWAS_frozen",
}


def save(frame: pd.DataFrame, out: Path, name: str) -> None:
    frame.to_csv(out / name, sep="\t", index=False)


def json_save(value: dict, path: Path) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def validate_boundary(frame: pd.DataFrame, registry: pd.DataFrame) -> None:
    """Require explicit dataset/sample membership in an authorized registry."""
    required = ["dataset_id", "sample_id"]
    if not set(required).issubset(frame) or not {*required, "authorized"}.issubset(registry):
        raise ValueError("Explicit dataset and sample membership/authorization are required.")
    if frame[required].isna().any().any() or registry[required].isna().any().any():
        raise ValueError("Missing dataset/sample identifiers are not allowed.")
    if not set(frame.dataset_id).issubset(ALLOWED):
        raise ValueError("Unapproved or mixed cohort input is excluded from this manuscript.")
    if registry.duplicated(required).any():
        raise ValueError("Registry contains ambiguous sample membership.")
    match = (
        frame[required]
        .drop_duplicates()
        .merge(registry, on=required, how="left", validate="one_to_one")
    )
    if not match.authorized.eq(True).all():
        raise ValueError("Sample absent from the authorized manuscript registry.")


def assert_hash(path: Path, expected: str) -> str:
    observed = file_sha256(path)
    if observed != expected:
        raise ValueError(f"Input checksum mismatch: {path.name}")
    return observed


def assert_reference(manifest: dict, expected: dict) -> None:
    for key, value in expected.items():
        if manifest.get(key) != value:
            raise ValueError(f"RFU backend/reference incompatibility: {key}")
    if manifest.get("failed_chunk_count") != 0 or manifest.get(
        "completed_chunk_count"
    ) != manifest.get("chunk_count"):
        raise ValueError("Assignment run is incomplete.")


def saved_rfu_labels(values: pd.Series) -> pd.Series:
    """Preserve source RFU labels; do not reinterpret IDs or shift numbering."""
    labels = values.astype("string")
    if labels.isna().any() or not labels.str.fullmatch(r"RFU[1-9]\d*").all():
        raise ValueError("Expected explicit saved one-based RFU labels.")
    return labels


def completed(out: Path, fingerprint: str) -> bool:
    path = out / "completion.json"
    if not path.exists():
        return False
    status = json.loads(path.read_text())
    if status.get("status") != "complete":
        return False
    if status.get("fingerprint") != fingerprint:
        raise ValueError(
            "Completed output has different scientific inputs/code/settings; use a new output directory."
        )
    for name, expected in status["outputs"].items():
        if not (out / name).is_file() or file_sha256(out / name) != expected:
            raise ValueError(f"Completed output missing or changed: {name}")
    return True


def geo_samples(path: Path) -> pd.DataFrame:
    records, record = [], None
    with gzip.open(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith("^SAMPLE = "):
                if record is not None:
                    records.append(record)
                record = {"accession": line.split(" = ", 1)[1], "files": [], "characteristics": []}
            elif record is not None and " = " in line:
                key, value = line.split(" = ", 1)
                if key == "!Sample_title":
                    record["title"] = value
                elif key == "!Sample_characteristics_ch1":
                    record["characteristics"].append(value)
                elif key.startswith("!Sample_supplementary_file") and value != "NONE":
                    record["files"].append(value.replace("ftp://", "https://"))
    if record is not None:
        records.append(record)
    return pd.DataFrame(records)


def reconcile_libraries(observations: pd.DataFrame, libraries: pd.DataFrame) -> pd.DataFrame:
    """Link processed prefixes to GEO by exact annotated donor/visit sets.

    Prefix numbers and GEO Batch numbers are not assumed equal. This is source
    metadata reconciliation, not independent demultiplexing or receptor matching.
    """
    result = []
    for prefix, rows in observations.groupby("library", dropna=False, sort=True):
        if pd.isna(prefix):
            raise ValueError("Missing processed-library prefix.")
        members = set(zip(rows.patient_rna, rows.state_rna, strict=True))
        for assay in ("RNA", "TCR"):
            candidates = []
            for geo_library, group in libraries[libraries.assay.eq(assay)].groupby("library"):
                if members == set(zip(group.patient_rna, group.state_rna, strict=True)):
                    candidates.append((geo_library, group.accession.unique()))
            if len(candidates) != 1 or len(candidates[0][1]) != 1:
                raise ValueError("Ambiguous or unsupported source patient/visit/library linkage.")
            result.append(
                {
                    "processed_library": prefix,
                    "geo_library": candidates[0][0],
                    "assay": assay,
                    "accession": candidates[0][1][0],
                    "linkage_basis": "exact source-annotated donor/visit membership; not receptor sequence",
                }
            )
    mapped = pd.DataFrame(result)
    if mapped.duplicated(["geo_library", "assay"]).any():
        raise ValueError("Library reconciliation is not bijective.")
    return mapped


def load_gse(root: Path, geo: Path, config: dict, out: Path) -> tuple[pd.DataFrame, dict]:
    raw = pd.read_csv(root / "GSE190905_TCR_data.csv.gz").rename(columns={"Unnamed: 0": "cell_id"})
    meta = pd.read_csv(root / "GSE190905_meta_data.csv.gz").rename(
        columns={"Unnamed: 0": "cell_id"}
    )
    rfu = pd.read_csv(root / "rfu/rfu_results_per_row.tsv.gz", sep="\t")
    if any(x.cell_id.isna().any() or x.cell_id.duplicated().any() for x in (raw, meta, rfu)):
        raise ValueError("Ambiguous source cell identifiers.")
    original_n = len(raw)
    frame = raw.merge(meta, on="cell_id", suffixes=("_tcr", "_rna"), validate="one_to_one")
    if len(frame) != original_n or set(rfu.cell_id) != set(raw.cell_id):
        raise ValueError("Source TCR/RNA/assignment coverage differs.")
    if (
        not frame.method_tcr.eq(frame.method_rna).all()
        or not frame.state_tcr.map({"before": "Pre", "after": "Post"}).eq(frame.state_rna).all()
    ):
        raise ValueError("Source treatment/timepoint conflict.")
    if (
        not frame.groupby("patient_tcr").patient_rna.nunique().eq(1).all()
        or not frame.groupby("patient_rna").patient_tcr.nunique().eq(1).all()
    ):
        raise ValueError("Source donor crosswalk is not one-to-one.")
    receptor_columns = {
        "TRB_1_cdr3": "cdr3aa",
        "TRB_1_cdr3_nt": "junction",
        "TRB_1_v_gene": "v_call",
        "TRB_1_j_gene": "j_call",
    }
    frame = frame.merge(
        rfu[["cell_id", *receptor_columns.values(), "rfu_label_nearest", "rfu_pass_threshold"]],
        on="cell_id",
        validate="one_to_one",
    )
    for original, cached in receptor_columns.items():
        if (
            frame[[original, cached]].isna().any().any()
            or not frame[original].eq(frame[cached]).all()
        ):
            raise ValueError(f"Saved assignment/source identity mismatch: {original}")
    if not frame.rfu_pass_threshold.isin([True, False]).all():
        raise ValueError("Assignment status must be explicit boolean.")
    donor_map = {p: f"D{i + 1:02}" for i, p in enumerate(sorted(meta.patient.unique()))}
    frame["donor"] = frame.patient_rna.map(donor_map)
    frame["time"] = frame.state_rna.map({"Pre": 0, "Post": 1})
    frame["sample_id"] = frame.donor + "_" + frame.state_rna
    frame["treatment"] = frame.method_rna
    frame["dataset_id"] = "GSE190905"
    frame["library"] = frame.cell_id.str.extract(r"^(b\d+)_", expand=False)
    # Physical libraries come from barcodes AND GEO; never from shared receptors.
    geo_rows = geo_samples(geo)
    library_rows = []
    for row in geo_rows.itertuples():
        match = re.fullmatch(r"(.+)_(RNA|TCR) \[Batch(\d+)\]", row.title)
        if not match:
            raise ValueError("Unexpected GEO library title; curate explicitly.")
        for member in match[1].split("/"):
            patient, visit = member.rsplit("_", 1)
            library_rows.append(
                {
                    "library": "b" + match[3],
                    "patient_rna": patient,
                    "state_rna": {"pre": "Pre", "post": "Post"}[visit],
                    "assay": match[2],
                    "accession": row.accession,
                }
            )
    libraries = pd.DataFrame(library_rows)
    library_meta = meta.rename(columns={"patient": "patient_rna", "state": "state_rna"})
    library_meta["library"] = library_meta.cell_id.str.extract(r"^(b\d+)_", expand=False)
    library_map = reconcile_libraries(library_meta, libraries)
    save(library_map, out, "gse190905_processed_geo_library_crosswalk.tsv")
    save(libraries, out, "gse190905_geo_library_membership.tsv")
    crosswalk = frame[["donor", "patient_rna", "patient_tcr"]].drop_duplicates()
    save(crosswalk, out, "gse190905_source_donor_crosswalk.tsv")
    # Primary TRB identity is not a paired alpha/beta clonotype.
    frame["clone"] = frame[["v_call", "j_call", "junction"]].agg("|".join, axis=1)
    frame["rfu"] = saved_rfu_labels(frame.rfu_label_nearest)
    frame["vj"] = frame.v_call + "|" + frame.j_call
    frame["length_bin"] = (frame.cdr3aa.str.len() // 5).astype(str)
    frame["compartment"] = "other"
    for label, states in config["gse190905_compartments"].items():
        frame.loc[frame.cluster.isin(states), "compartment"] = label
    columns = ["rfu", "rfu_pass_threshold", "v_call", "j_call", "vj", "length_bin", "cdr3aa"]
    if frame.groupby("clone")[columns].nunique(dropna=False).gt(1).any().any():
        raise ValueError("A receptor identity does not have a fixed group/feature mapping.")
    manifest = (
        frame.groupby(["dataset_id", "sample_id", "donor", "time", "treatment", "library"])
        .agg(
            tcr_cells=("cell_id", "size"),
            unique_primary_trb=("clone", "nunique"),
            threshold_cells=("rfu_pass_threshold", "sum"),
        )
        .reset_index()
    )
    meta["donor"] = meta.patient.map(donor_map)
    meta["sample_id"] = meta.donor + "_" + meta.state
    manifest["rna_metadata_cells"] = manifest.sample_id.map(meta.groupby("sample_id").size())
    manifest["tcr_coverage"] = manifest.tcr_cells / manifest.rna_metadata_cells
    manifest["threshold_coverage_of_tcr"] = manifest.threshold_cells / manifest.tcr_cells
    manifest["authorized"] = True
    validate_boundary(frame, manifest[["dataset_id", "sample_id", "authorized"]])
    design = validate_longitudinal_design(
        frame, sample_key="sample_id", donor_key="donor", time_key="time", condition_key="treatment"
    )
    save(manifest, out, "gse190905_paired_manifest.tsv")
    save(design.qc_table, out, "gse190905_design_qc.tsv")
    state = (
        frame.groupby(["donor", "sample_id", "time", "treatment", "cluster"])
        .size()
        .rename("tcr_cells")
        .reset_index()
    )
    state["fraction_of_tcr_cells"] = state.tcr_cells / state.groupby(
        "sample_id"
    ).tcr_cells.transform("sum")
    save(state, out, "gse190905_state_composition.tsv")
    return frame, {
        "tcr_cells": len(frame),
        "rna_metadata_cells": len(meta),
        "donors": frame.donor.nunique(),
        "samples": frame.sample_id.nunique(),
        "physical_libraries": frame.library.nunique(),
        "source_treatment_donors": frame.drop_duplicates("donor")
        .treatment.value_counts()
        .to_dict(),
        "library_linkage": "GEO library titles reconciled to barcode-prefix pools by exact RNA donor/time membership; processed prefix numbers differ from GEO Batch numbers",
        "final_paper_discrepancy": "cached six-donor release, not final seven-donor cohort; GEO metastatic wording conflicts with paper early-stage wording",
        "cell_state_denominator": "TCR-bearing cells only; RNA-only cells lack source cluster labels in metadata table",
    }


def pair_counts(frame: pd.DataFrame, weighting: str) -> tuple[pd.Series, pd.Series]:
    result = []
    for time in (0, 1):
        counts = frame.loc[frame.time.eq(time)].groupby("clone").size()
        result.append(counts if weighting == "cell" else counts.clip(upper=1))
    return tuple(result)


def paired_measurements(frame: pd.DataFrame, config: dict, out: Path) -> dict:
    table, group_rows, persistence, support = [], [], [], []
    subsets = {
        "all_T": frame,
        **{c: frame[frame.compartment.eq(c)] for c in ("CD4", "CD8")},
        **{"state:" + s: frame[frame.cluster.eq(s)] for s in sorted(frame.cluster.unique())},
    }
    mapping = frame.drop_duplicates("clone").set_index("clone")
    for name, subset in subsets.items():
        for donor in sorted(frame.donor.unique()):
            person = subset[subset.donor.eq(donor)]
            for policy in ("threshold", "nearest"):
                kept = person[person.rfu_pass_threshold] if policy == "threshold" else person
                sizes = kept.groupby("time").size().reindex([0, 1], fill_value=0)
                eligible = bool(sizes.ge(config["min_cells_per_visit"]).all())
                support.append(
                    {
                        "donor": donor,
                        "subset": name,
                        "policy": policy,
                        "before_cells": sizes[0],
                        "after_cells": sizes[1],
                        "status": "analyzed"
                        if eligible
                        else "insufficient_observed_cells_or_missing_visit",
                    }
                )
                if not eligible:
                    continue
                total = person.groupby("time").size()
                for weighting in config["weights"]:
                    before, after = pair_counts(kept, weighting)
                    base = {
                        "donor": donor,
                        "treatment": kept.treatment.iloc[0],
                        "subset": name,
                        "policy": policy,
                        "weighting": weighting,
                        "cell_coverage_before": sizes[0] / total[0],
                        "cell_coverage_after": sizes[1] / total[1],
                    }
                    for name_col, representation in (
                        ("rfu", "RFU"),
                        ("v_call", "TRBV"),
                        ("vj", "TRBV_TRBJ"),
                    ):
                        change = multiscale_repertoire_change(before, after, mapping[name_col])
                        table.append({**base, "grouping": representation, **change.summary})
                        if representation == "RFU":
                            g = change.group_changes.reset_index().rename(
                                columns={"group": "rfu_label"}
                            )
                            group_rows.append(g.assign(**base))
                    if weighting == "cell":
                        for rfu, members in kept.groupby("rfu", sort=True):
                            a, b = pair_counts(members, "cell")
                            for detection in config["rfu_detection_cell_counts"]:
                                persistence.append(
                                    {
                                        **base,
                                        "rfu_label": rfu,
                                        "detection_min_cells": detection,
                                        "before_cells": int(a.sum()),
                                        "after_cells": int(b.sum()),
                                        "observed_before": bool(a.sum() >= detection),
                                        "observed_after": bool(b.sum() >= detection),
                                        "before_unique_clones": len(a),
                                        "after_unique_clones": len(b),
                                        "observed_shared_clones": len(
                                            a.index.intersection(b.index)
                                        ),
                                        "dominant_clone_fraction_before": a.max() / a.sum()
                                        if len(a)
                                        else np.nan,
                                        "dominant_clone_fraction_after": b.max() / b.sum()
                                        if len(b)
                                        else np.nan,
                                        "one_cell_fraction_before": 1 / sizes[0],
                                        "one_cell_fraction_after": 1 / sizes[1],
                                    }
                                )
    measurements = pd.DataFrame(table)
    save(measurements, out, "gse190905_multiscale_pairs.tsv")
    save(pd.concat(group_rows, ignore_index=True), out, "gse190905_rfu_changes.tsv")
    save(pd.DataFrame(persistence), out, "gse190905_rfu_observed_persistence.tsv")
    save(pd.DataFrame(support), out, "gse190905_pair_support.tsv")
    # Similarity summaries give each donor one within and one averaged-between value.
    selected = frame[frame.rfu_pass_threshold]
    matrix = pd.crosstab(selected.sample_id, selected.rfu)
    metadata = selected[["sample_id", "donor", "time"]].drop_duplicates().set_index("sample_id")
    pairs = longitudinal_similarity(
        matrix, metadata=metadata, donor_key="donor", time_key="time", metric="cosine"
    )
    save(pairs, out, "gse190905_rfu_similarity_pairs.tsv")
    similarities = []
    for donor in sorted(selected.donor.unique()):
        for relation in (True, False):
            part = pairs[
                pairs.same_donor.eq(relation) & (pairs.donor_a.eq(donor) | pairs.donor_b.eq(donor))
            ]
            similarities.append(
                {
                    "donor": donor,
                    "relation": "within" if relation else "between_mean",
                    "cosine": part.value.mean(),
                    "n_pairs_contributing": len(part),
                }
            )
    save(pd.DataFrame(similarities), out, "gse190905_rfu_similarity_donor_summary.tsv")
    primary = measurements[
        (measurements.subset == "all_T")
        & (measurements.policy == "threshold")
        & (measurements.weighting == "cell")
    ]
    summary = {
        label: {
            "n_donors": len(x),
            "d_clone_median": float(x.d_clone.median()),
            "d_group_median": float(x.d_group.median()),
            "cancellation_median": float(x.aggregation_cancellation.median()),
        }
        for label, x in primary.groupby("grouping")
    }
    return summary


def sensitivities(frame: pd.DataFrame, config: dict, out: Path) -> None:
    selected = frame[frame.rfu_pass_threshold].copy()
    mapping = selected.drop_duplicates("clone").set_index("clone")
    strata = mapping[["v_call", "j_call", "length_bin"]]
    variations = (
        strata.assign(rfu=mapping.rfu).groupby(list(strata.columns)).rfu.transform("nunique")
    )
    fixed_random = [
        permute_fixed_groups(mapping.rfu, strata=strata, random_state=config["random_seed"] + i)
        for i in range(config["random_group_replicates"])
    ]
    controls, sampling, dominance = [], [], []
    for donor_index, (donor, person) in enumerate(selected.groupby("donor", sort=True)):
        before, after = pair_counts(person, "cell")
        if not len(before) or not len(after):
            continue
        for i, grouping in enumerate(fixed_random):
            x = multiscale_repertoire_change(before, after, grouping).summary
            controls.append(
                {
                    "donor": donor,
                    "replicate": i,
                    "seed": config["random_seed"] + i,
                    "global_group_sizes_and_feature_strata_preserved": True,
                    "fraction_receptors_in_exchangeable_strata": float(variations.gt(1).mean()),
                    "fraction_labels_changed": float(grouping.ne(mapping.rfu).mean()),
                    **x,
                }
            )
        cells_before = person.loc[person.time.eq(0), "clone"].to_numpy()
        cells_after = person.loc[person.time.eq(1), "clone"].to_numpy()
        n = min(config["depth_cap_cells"], len(cells_before), len(cells_after))
        for i in range(config["depth_replicates"]):
            seed = config["random_seed"] + 10000 * (donor_index + 1) + i
            rng = np.random.default_rng(seed)
            # Actual cells without replacement: not multinomial simulation.
            a = pd.Series(rng.choice(cells_before, n, replace=False)).value_counts()
            b = pd.Series(rng.choice(cells_after, n, replace=False)).value_counts()
            change = multiscale_repertoire_change(a, b, mapping.rfu).summary
            sampling.append(
                {
                    "donor": donor,
                    "replicate": i,
                    "seed": seed,
                    "cells_per_visit": n,
                    "scheme": "observed_cell_subsampling_without_replacement",
                    **change,
                }
            )
        removed = {before.idxmax(), after.idxmax()}
        a, b = (
            before.drop(list(removed), errors="ignore"),
            after.drop(list(removed), errors="ignore"),
        )
        change = multiscale_repertoire_change(a, b, mapping.rfu).summary
        dominance.append(
            {
                "donor": donor,
                "scheme": "remove_union_of_each_visit_dominant_primary_TRB",
                "removed_clones": len(removed),
                "retained_cell_fraction_before": a.sum() / before.sum(),
                "retained_cell_fraction_after": b.sum() / after.sum(),
                **change,
            }
        )
    save(pd.DataFrame(controls), out, "gse190905_fixed_group_controls.tsv")
    save(pd.DataFrame(sampling), out, "gse190905_empirical_depth_sensitivity.tsv")
    save(pd.DataFrame(dominance), out, "gse190905_dominant_clone_sensitivity.tsv")


def wells_context(root: Path, config: dict, out: Path) -> dict:
    meta = pd.read_csv(root / "full_run/obs_metadata.tsv.gz", sep="\t")
    rfu = pd.read_csv(root / "full_run/rfu_results_per_row.tsv.gz", sep="\t")
    frame = rfu.merge(meta, on="cell_id", validate="one_to_one")
    if len(frame) != len(rfu) or frame.donor_id.isna().any():
        raise ValueError("Wells source metadata/assignment coverage is incomplete.")
    # Include atlas donors with no eligible TRB so coverage cannot hide them.
    donor_map = {x: f"W{i + 1:02}" for i, x in enumerate(sorted(meta.donor_id.unique()))}
    frame["donor"] = frame.donor_id.map(donor_map)
    frame["dataset_id"] = "Wells_public"
    frame["sample_id"] = frame.library_id
    registry = frame[["dataset_id", "sample_id"]].drop_duplicates().assign(authorized=True)
    validate_boundary(frame, registry)
    save(
        pd.DataFrame({"source_donor": list(donor_map), "donor": list(donor_map.values())}),
        out,
        "wells_source_donor_crosswalk.tsv",
    )
    coverage = (
        frame.groupby(["donor", "tissue"])
        .agg(primary_trb_cells=("cell_id", "size"), threshold_cells=("rfu_pass_threshold", "sum"))
        .reset_index()
    )
    meta["donor"] = meta.donor_id.map(donor_map)
    all_cells = meta.groupby(["donor", "tissue"]).size().rename("atlas_cells").reset_index()
    coverage = all_cells.merge(
        coverage, on=["donor", "tissue"], how="left", validate="one_to_one"
    ).fillna(0)
    coverage["primary_trb_coverage"] = coverage.primary_trb_cells / coverage.atlas_cells
    coverage["threshold_fraction_of_primary_trb"] = (
        coverage.threshold_cells / coverage.primary_trb_cells.replace(0, np.nan)
    )
    save(coverage, out, "wells_donor_tissue_coverage.tsv")
    selected = frame[frame.rfu_pass_threshold].copy()
    selected["rfu_label"] = saved_rfu_labels(selected.rfu_label_nearest)
    selected["observed_receptor"] = selected.cdr3aa + "|" + selected.v_call.fillna("V_UNKNOWN")
    group = ["donor", "tissue", "rfu_label", "cell_type"]
    context = (
        selected.groupby(group)
        .agg(cells=("cell_id", "size"), unique_observed_receptors=("observed_receptor", "nunique"))
        .reset_index()
    )
    context["fraction_within_donor_tissue_rfu"] = context.cells / context.groupby(
        group[:-1]
    ).cells.transform("sum")
    save(context, out, "wells_rfu_donor_tissue_celltype.tsv")
    support = (
        selected.groupby("rfu_label")
        .agg(
            cells=("cell_id", "size"),
            donors=("donor", "nunique"),
            tissues=("tissue", "nunique"),
            observed_receptors=("observed_receptor", "nunique"),
        )
        .reset_index()
    )
    support["public_support_rule"] = support.cells.ge(
        config["wells_min_cells"]
    ) & support.donors.ge(config["wells_min_donors"])
    save(support, out, "wells_rfu_support.tsv")
    dominance = (
        selected.groupby(["donor", "tissue", "rfu_label", "observed_receptor"])
        .size()
        .rename("cells")
        .reset_index()
    )
    dominance = (
        dominance.groupby(["donor", "tissue", "rfu_label"])
        .agg(
            cells=("cells", "sum"),
            unique_receptors=("cells", "size"),
            dominant_clone_cells=("cells", "max"),
        )
        .reset_index()
    )
    dominance["dominant_observed_receptor_fraction"] = (
        dominance.dominant_clone_cells / dominance.cells
    )
    save(dominance, out, "wells_rfu_clonal_dominance.tsv")
    return {
        "atlas_cells": len(meta),
        "atlas_donors": meta.donor_id.nunique(),
        "primary_trb_cells": len(frame),
        "threshold_cells": len(selected),
        "donors": frame.donor.nunique(),
        "tissues": frame.tissue.nunique(),
        "threshold_rfus": len(support),
        "rfus_meeting_public_support_rule": int(support.public_support_rule.sum()),
        "assignment_reuse": "standard backend; no dependence on historical map-aware candidate subset",
        "context_label": "source cell_type (cell_state is mostly 'na'); cell counts are not bulk template/read counts",
        "clone_limit": "observed CDR3 amino acid + V, not paired-chain/nucleotide lineage",
        "old_map_aware_parity": "not established; old cache absent and not used",
    }


def external_inventory(geo: Path, out: Path) -> dict:
    samples = geo_samples(geo)
    rows = []
    for row in samples.itertuples():
        match = re.fullmatch(r"HyPR-HN (\d+)_(\d+)(?:_(PBMC))?(?: tumor)?, (GEX|TCR)", row.title)
        for url in row.files or [""]:
            rows.append(
                {
                    "accession": row.accession,
                    "title": row.title,
                    "longitudinal_hypr_sample": bool(match),
                    "donor": match[1] if match else "",
                    "visit": match[2] if match else "",
                    "compartment": ("blood" if match[3] else "tumor") if match else "",
                    "assay": match[4] if match else "other",
                    "url": url,
                    "characteristics": "; ".join(row.characteristics),
                    "downloaded_for_endpoints": False,
                }
            )
    files = pd.DataFrame(rows)
    save(files, out, "gse280982_processed_file_inventory.tsv")
    subjects = files[files.longitudinal_hypr_sample].drop_duplicates(
        ["donor", "visit", "compartment", "assay"]
    )
    pairs = subjects.pivot_table(
        index=["donor", "visit", "compartment"],
        columns="assay",
        values="accession",
        aggfunc="first",
    ).reset_index()
    pairs["paired_processed_GEX_TCR_available"] = pairs[["GEX", "TCR"]].notna().all(axis=1)
    save(pairs, out, "gse280982_public_pair_availability.tsv")
    return {
        "geo_samples_all_assays": len(samples),
        "hypr_paired_processed_visits": int(pairs.paired_processed_GEX_TCR_available.sum()),
        "paired_tumor_visits": int(
            (pairs.paired_processed_GEX_TCR_available & pairs.compartment.eq("tumor")).sum()
        ),
        "tumor_TCR_donors": subjects[
            (subjects.assay == "TCR") & (subjects.compartment == "tumor")
        ].donor.nunique(),
        "status": "processed contigs and GEX matrices listed; no new assignment or endpoint analysis this checkpoint",
        "prior_access": "GEO metadata and paper inspected before external endpoints; not claimed untouched",
    }


def figures(out: Path, fingerprint: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    set_style()
    data = pd.read_csv(out / "gse190905_multiscale_pairs.tsv", sep="\t")
    primary = data[
        (data.subset == "all_T") & (data.policy == "threshold") & (data.weighting == "cell")
    ]
    fig, axes = plt.subplots(1, 3, figsize=(9, 2.8), layout="constrained")
    colors = plt.get_cmap("tab10")
    metrics = ["d_clone", "RFU", "TRBV", "TRBV_TRBJ"]
    for i, (donor, x) in enumerate(primary.groupby("donor", sort=True)):
        by_group = x.set_index("grouping")
        axes[0].plot(
            range(4),
            [x.d_clone.iloc[0], *[by_group.loc[k, "d_group"] for k in metrics[1:]]],
            "o-",
            color=colors(i),
            label=donor,
        )
    axes[0].set(
        xticks=range(4),
        xticklabels=["TRB clone", "RFU", "TRBV", "V/J"],
        ylabel="Pre/post total variation",
        ylim=(0, 1.05),
        title="Six paired donors; assigned cells",
    )
    coverage = pd.read_csv(out / "gse190905_paired_manifest.tsv", sep="\t")
    for i, (donor, x) in enumerate(coverage.groupby("donor", sort=True)):
        x = x.sort_values("time")
        axes[1].plot(x.time, x.threshold_coverage_of_tcr, "o-", color=colors(i), label=donor)
    axes[1].set(
        xticks=[0, 1],
        xticklabels=["Pre", "Post"],
        ylabel="Threshold coverage of TCR cells",
        ylim=(0, 1.05),
        title="Coverage remains a separate endpoint",
    )
    sims = pd.read_csv(out / "gse190905_rfu_similarity_donor_summary.tsv", sep="\t")
    for i, (_donor, x) in enumerate(sims.groupby("donor", sort=True)):
        values = x.set_index("relation").cosine
        axes[2].plot([0, 1], values.reindex(["within", "between_mean"]), "o-", color=colors(i))
    axes[2].set(
        xticks=[0, 1],
        xticklabels=["Within", "Between mean"],
        ylabel="RFU cosine similarity",
        ylim=(0, 1.05),
        title="One paired summary per donor",
    )
    axes[0].legend(ncol=3, frameon=False, loc="lower left")
    for label, ax in zip("ABC", axes, strict=True):
        clean_axis(ax)
        panel_label(ax, label)
    save_figure(fig, out, "gse190905_longitudinal")
    plt.close(fig)
    depth = pd.read_csv(out / "gse190905_empirical_depth_sensitivity.tsv", sep="\t")
    control = pd.read_csv(out / "gse190905_fixed_group_controls.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(7, 2.8), layout="constrained")
    donors = sorted(primary.donor.unique())
    for i, donor in enumerate(donors):
        observed = primary[(primary.donor == donor) & (primary.grouping == "RFU")].iloc[0]
        for ax, source, metric in [
            (axes[0], depth, "d_group"),
            (axes[1], control, "aggregation_cancellation"),
        ]:
            values = source.loc[source.donor.eq(donor), metric]
            q = values.quantile([0.025, 0.5, 0.975])
            ax.errorbar(
                i,
                q.loc[0.5],
                yerr=[[q.loc[0.5] - q.loc[0.025]], [q.loc[0.975] - q.loc[0.5]]],
                fmt="o",
                color="#777777",
                capsize=3,
            )
            ax.scatter(i, observed[metric], marker="D", color="#0072B2", zorder=3)
    axes[0].set(title="Observed-cell depth sensitivity", ylabel="RFU total variation")
    axes[1].set(title="Fixed feature/size-matched groupings", ylabel="Aggregation cancellation")
    for label, ax in zip("AB", axes, strict=True):
        ax.set_xticks(range(len(donors)), donors)
        clean_axis(ax)
        panel_label(ax, label)
    fig.suptitle(
        "Blue: observed RFU; gray: median and 95% replicate range (not donor CI)", fontsize=8
    )
    save_figure(fig, out, "gse190905_multiscale_sensitivity")
    plt.close(fig)
    panels = [
        (
            "gse190905_longitudinal",
            "A",
            "gse190905_multiscale_pairs.tsv",
            "six donors",
            "threshold-assigned TCR cells; same cells for all groupings",
            "primary TRB clone != paired alpha/beta clone",
        ),
        (
            "gse190905_longitudinal",
            "B",
            "gse190905_paired_manifest.tsv",
            "six donors",
            "TCR cells at each visit",
            "no absolute depletion inference",
        ),
        (
            "gse190905_longitudinal",
            "C",
            "gse190905_rfu_similarity_donor_summary.tsv",
            "six donors",
            "one within and one averaged-between value per donor",
            "shared between-donor comparisons are dependent",
        ),
        (
            "gse190905_multiscale_sensitivity",
            "A",
            "gse190905_empirical_depth_sensitivity.tsv",
            "six donors",
            "actual assigned cells without replacement",
            "replicate intervals are not participant uncertainty",
        ),
        (
            "gse190905_multiscale_sensitivity",
            "B",
            "gse190905_fixed_group_controls.tsv",
            "six donors",
            "global receptor-feature strata and group sizes",
            "descriptive control; no exchangeability or calibrated p-value",
        ),
    ]
    source = []
    for figure, panel, table, unit, denominator, limitation in panels:
        tables = [table]
        if figure == "gse190905_multiscale_sensitivity":
            tables.append("gse190905_multiscale_pairs.tsv")
        source.append(
            {
                "figure": figure,
                "panel": panel,
                "source_table": ";".join(tables),
                "source_sha256": ";".join(file_sha256(out / name) for name in tables),
                "input_config_code_fingerprint": fingerprint,
                "generating_command": "python -m manuscript.scripts.radiation_methods (full command in provenance.json)",
                "biological_unit": unit,
                "denominator": denominator,
                "interpretation_limit": limitation,
                "completion_status": "complete",
            }
        )
    save(pd.DataFrame(source), out, "figure_source_index.tsv")


def run(args: argparse.Namespace) -> None:
    config = json.loads(args.config.read_text())
    if not set(config["allowed_datasets"]).issubset(ALLOWED):
        raise ValueError("Configuration attempts to admit an unapproved cohort.")
    args.out.mkdir(parents=True, exist_ok=True)
    pinned = {
        "GSE190905_TCR_data.csv.gz": args.gse / "GSE190905_TCR_data.csv.gz",
        "GSE190905_meta_data.csv.gz": args.gse / "GSE190905_meta_data.csv.gz",
        "GSE190905_per_row": args.gse / "rfu/rfu_results_per_row.tsv.gz",
        "Wells_per_row": args.wells / "full_run/rfu_results_per_row.tsv.gz",
    }
    hashes = {key: assert_hash(path, config["source_sha256"][key]) for key, path in pinned.items()}
    paths = {
        **pinned,
        "config": args.config,
        "Wells_metadata": args.wells / "full_run/obs_metadata.tsv.gz",
        "Wells_run": args.wells / "full_run/run_manifest.json",
        "GSE190905_run": args.gse / "rfu/run_manifest.json",
        "Wells_evidence": args.wells / "evidence_manifest_full_rfu.json",
        "GSE190905_evidence": args.gse / "evidence_manifest_rfu.json",
        "GSE190905_GEO": args.sources / "GSE190905_family.soft.gz",
        "GSE280982_GEO": args.sources / "GSE280982_family.soft.gz",
    }
    for path in [
        Path(__file__),
        Path(__file__).with_name("_figure_common.py"),
        REPO / "src/scrfu/multiscale.py",
        REPO / "src/scrfu/longitudinal.py",
        REPO / "src/scrfu/completed_run.py",
    ]:
        paths[str(path.relative_to(REPO))] = path
    hashes.update({key: file_sha256(path) for key, path in paths.items()})
    versions = {
        name: importlib.metadata.version(name)
        for name in ["numpy", "pandas", "matplotlib", "scrfu", "anndata"]
    }
    fingerprint = hashlib.sha256(
        json.dumps(
            {"hashes": hashes, "packages": versions, "python": sys.version}, sort_keys=True
        ).encode()
    ).hexdigest()
    if completed(args.out, fingerprint):
        print("Completed manuscript stage verified; reused without recomputation.")
        return
    json_save(
        {
            "status": "running",
            "fingerprint": fingerprint,
            "scope": "public development and context; RP1-14 blocked on asset location",
        },
        args.out / "completion.json",
    )
    evidence = []
    for dataset, root, evidence_manifest in [
        ("GSE190905", args.gse / "rfu", args.gse / "evidence_manifest_rfu.json"),
        ("Wells_public", args.wells / "full_run", args.wells / "evidence_manifest_full_rfu.json"),
    ]:
        manifest = json.loads((root / "run_manifest.json").read_text())
        assert_reference(manifest, config["reference"])
        validation = validate_completed_rfu_run(
            root, evidence_manifest=evidence_manifest, expected_provenance=config["reference"]
        )
        json_save(validation, args.out / f"{dataset}_reuse_validation.json")
        evidence.append(
            {
                "dataset": dataset,
                "current_path": str(root),
                "stage": "completed primary-TRB assignment",
                "values": "one primary TRB per cell; no bulk read counts",
                "fields": "donor/sample/visit/treatment (GSE190905) or donor/library/tissue/cell_type (Wells)",
                "backend": manifest["backend_mode"],
                "reference_sha256": manifest["km5000_rdata_sha256"],
                "threshold": manifest["rfu_threshold"],
                "generating_script": manifest["workflow"],
                "integrity": "recorded assignment hashes and completed-run validator",
                "access_status": "public-derived; original participant labels not committed",
                "disposition": "direct reuse",
            }
        )
    evidence.append(
        {
            "dataset": "RP1-14_published_authorized",
            "current_path": "not located in bounded historical-path/inventory checks",
            "stage": "historical assignment known from user; matrix unavailable",
            "access_status": "published/authorized scope still requires recovered sample registry; not redistributable by assumption",
            "disposition": "blocked; do not substitute pooled RP/RERF transforms",
        }
    )
    save(pd.DataFrame(evidence), args.out, "local_asset_inventory.tsv")
    frame, gse = load_gse(args.gse, paths["GSE190905_GEO"], config, args.out)
    print("GSE190905 metadata, source identity and library linkage validated.", flush=True)
    measures = paired_measurements(frame, config, args.out)
    sensitivities(frame, config, args.out)
    print("Paired multiscale measurements and fixed sensitivity analyses complete.", flush=True)
    wells = wells_context(args.wells, config, args.out)
    external = external_inventory(paths["GSE280982_GEO"], args.out)
    figures(args.out, fingerprint)
    counts = {
        "GSE190905": gse,
        "primary_multiscale": measures,
        "Wells_public": wells,
        "GSE280982": external,
        "RP1-14": {
            "status": "blocked_missing_assets",
            "analysis_not_reconstructed_from_unpublished_cohort": True,
        },
    }
    json_save(counts, args.out / "evidence_counts.json")
    provenance = {
        "command": [sys.executable, "-m", "manuscript.scripts.radiation_methods", *sys.argv[1:]],
        "python": sys.executable,
        "scrfu_import": scrfu.__file__,
        "git_sha": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
        ).strip(),
        "git_dirty_at_execution": bool(
            subprocess.check_output(["git", "status", "--porcelain"], cwd=REPO, text=True).strip()
        ),
        "packages": versions,
        "inputs_and_code": {
            key: {"path": str(path), "sha256": hashes[key]} for key, path in paths.items()
        },
        "config": config,
        "fingerprint": fingerprint,
        "assignment_campaign": "none; reused completed runs",
        "regulatory": "frozen d62f1c1; no refit/download/reanalysis and no unverified RFU-label crosslink",
    }
    json_save(provenance, args.out / "provenance.json")
    names = sorted(
        p.name for p in args.out.iterdir() if p.is_file() and p.name != "completion.json"
    )
    json_save(
        {
            "status": "complete",
            "scope": "GSE190905 development, Wells context, external-file audit",
            "blocked_stages": [
                "RP1-14: missing local matrix/metadata",
                "GSE280982 endpoint analysis: not run this checkpoint",
            ],
            "fingerprint": fingerprint,
            "outputs": {name: file_sha256(args.out / name) for name in names},
        },
        args.out / "completion.json",
    )
    print(json.dumps(counts, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=REPO / "manuscript/config/radiation_methods_v1.json"
    )
    parser.add_argument("--gse", type=Path, required=True)
    parser.add_argument("--wells", type=Path, required=True)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
