"""Bounded Wells context preparation and render-only reusable figure exports.

No assignment, embedding fitting, or prediction model is called. Preparation
does one explicitly specified, exploratory donor-paired expression contrast.
Render iterations consume its frozen tables and never repeat that contrast.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import shlex
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import stats

from scrfu.io import read_h5ad_obs

from .rp1_14_reusable_figures import REPO, json_save, require, safe_child, sha256

SPEC = REPO / "manuscript/figures/singlecell_context_spec.json"
AA = "ACDEFGHIKLMNPQRSTVWY"
DOC = REPO / "manuscript/figures/singlecell_context_captions.md"


def read_table(path):
    return pd.read_csv(path, sep="\t", float_precision="round_trip", keep_default_na=False)


def save_table(frame, out, name):
    frame.to_csv(out / "tables" / f"{name}.tsv", sep="\t", index=False, float_format="%.17g")


def stable_order(ids, salt):
    return sorted(ids, key=lambda x: hashlib.sha256(f"{salt}|{x}".encode()).digest())


def read_column(group, name):
    obj = group[name]
    if isinstance(obj, h5py.Group):
        codes = obj["codes"][:]
        require((codes >= 0).all(), f"Missing categorical feature: {name}")
        categories = obj["categories"].asstr()[:]
        return categories[codes]
    return obj.asstr()[:] if obj.dtype.kind in {"O", "S", "U"} else obj[:]


def verify_sources(spec):
    for key, record in spec["sources"].items():
        require(sha256(Path(record["path"])) == record["sha256"], f"Changed frozen source: {key}")


def load_cells(spec):
    atlas = Path(spec["sources"]["atlas"]["path"])
    obs = read_h5ad_obs(atlas, columns=["cell_type", "tissue", "donor_id", "library_id"])
    require(obs.index.is_unique, "Duplicate source cell IDs")
    obs = obs.rename_axis("cell_id").reset_index()
    obs["source_row"] = np.arange(len(obs))
    aliases = {x: f"W{i + 1:02}" for i, x in enumerate(sorted(obs.donor_id.unique()))}
    obs["donor"] = obs.donor_id.map(aliases)
    assignment = pd.read_csv(
        spec["sources"]["assignments"]["path"],
        sep="\t",
        usecols=[
            "cell_id",
            "cdr3aa",
            "v_call",
            "rfu_label_nearest",
            "rfu_score",
            "rfu_pass_threshold",
        ],
    )
    require(assignment.cell_id.is_unique, "Multiple primary assignments per cell")
    require(assignment.cell_id.isin(obs.cell_id).all(), "Assignment cell absent from source atlas")
    require(
        assignment.rfu_pass_threshold.eq(assignment.rfu_score.ge(0.6)).all(), "Threshold flag drift"
    )
    cells = obs.merge(assignment, on="cell_id", how="left", validate="one_to_one", sort=False)
    cells["assignment_status"] = np.select(
        [cells.rfu_pass_threshold.eq(True), cells.rfu_score.notna()],
        ["qualified", "below_threshold"],
        default="no_primary_TRB",
    )
    cells["rfu"] = cells.rfu_label_nearest.where(cells.assignment_status.eq("qualified"), "")
    with h5py.File(atlas, "r") as f:
        xy = f[spec["embedding"]["key"]][:]
        require(xy.shape == (len(cells), 2) and np.isfinite(xy).all(), "Invalid frozen embedding")
        cells[["umap_1", "umap_2"]] = xy.astype(np.float64)
    return cells


def rank_rfus(qualified):
    counts = qualified.groupby(["rfu", "donor"], observed=True).size().unstack(fill_value=0)
    result = pd.DataFrame(
        {
            "cells": counts.sum(axis=1),
            "donors": counts.gt(0).sum(axis=1),
            "donors_ge5": counts.ge(5).sum(axis=1),
            "donors_ge10": counts.ge(10).sum(axis=1),
            "max_donor_cells": counts.max(axis=1),
        }
    ).reset_index()
    result = result.sort_values(
        ["donors_ge10", "donors_ge5", "cells", "rfu"], ascending=[False, False, False, True]
    ).reset_index(drop=True)
    result["selection_rank"] = np.arange(1, len(result) + 1)
    return result


def logo_tables(qualified, labels):
    sequences, frequencies, summaries = [], [], []
    for label in labels:
        cells = qualified[qualified.rfu.eq(label)]
        seq = (
            cells.groupby("cdr3aa", observed=True)
            .agg(
                cells=("cell_id", "size"),
                donors=("donor", "nunique"),
                v_calls=("v_call", lambda x: ";".join(sorted(set(x)))),
            )
            .reset_index()
        )
        require(seq.cdr3aa.map(lambda x: set(x) <= set(AA)).all(), "Unexpected logo amino acid")
        seq["length"] = seq.cdr3aa.str.len()
        length_counts = seq.groupby("length").size()
        length = int(length_counts[length_counts.eq(length_counts.max())].index.min())
        seq["rfu"] = label
        seq["included_in_logo"] = seq.length.eq(length)
        used = seq[seq.included_in_logo]
        for pos in range(length):
            counts = used.cdr3aa.str[pos].value_counts()
            for aa in AA:
                n = int(counts.get(aa, 0))
                frequencies.append(
                    {
                        "rfu": label,
                        "position": pos + 1,
                        "amino_acid": aa,
                        "sequence_count": n,
                        "denominator": len(used),
                        "frequency": n / len(used),
                    }
                )
        summaries.append(
            {
                "rfu": label,
                "cells": len(cells),
                "donors": cells.donor.nunique(),
                "unique_aa": len(seq),
                "modal_length": length,
                "logo_unique_aa": len(used),
                "logo_cells": int(used.cells.sum()),
                "max_donor_cells": int(cells.donor.value_counts().max()),
            }
        )
        sequences.append(seq)
    return (
        pd.concat(sequences, ignore_index=True),
        pd.DataFrame(frequencies),
        pd.DataFrame(summaries),
    )


def match_cells(qualified, contrast, seed):
    """Exact strata, balanced without replacement; only metadata enters selection."""
    data = qualified.copy()
    data["group"] = np.where(data.rfu.eq(contrast["case_rfu"]), "RFU_positive", "background")
    selected, strata = [], []
    keys = contrast["match_columns"]
    require(data[keys].notna().all().all(), "Missing matching metadata")
    for key, part in data.groupby(keys, observed=True, sort=True):
        pos, neg = part[part.group.eq("RFU_positive")], part[part.group.eq("background")]
        if pos.empty:
            continue
        n = min(len(pos), len(neg))
        identifier = f"S{len(strata) + 1:04}"
        strata.append(
            dict(zip(keys, key, strict=True))
            | {
                "stratum": identifier,
                "case_available": len(pos),
                "background_available": len(neg),
                "matched_per_group": n,
                "donor": part.donor.iloc[0],
            }
        )
        for group, source in [("RFU_positive", pos), ("background", neg)]:
            ids = stable_order(source.cell_id.tolist(), f"{seed}|{group}")[:n]
            selected.append(
                source.set_index("cell_id").loc[ids].reset_index().assign(stratum=identifier)
            )
    matched = pd.concat(selected, ignore_index=True)
    counts = matched.groupby(["donor", "group"], observed=True).size().unstack(fill_value=0)
    counts = counts.reindex(columns=["RFU_positive", "background"], fill_value=0)
    require(counts.RFU_positive.eq(counts.background).all(), "Unequal matched groups")
    retained = counts.index[counts.RFU_positive.ge(contrast["minimum_cells_per_group_per_donor"])]
    matched["included_in_de"] = matched.donor.isin(retained)
    strata = pd.DataFrame(strata)
    strata["donor_retained"] = strata.donor.isin(retained)
    require(matched.cell_id.is_unique, "A matched cell was reused")
    require(len(retained) >= 6, "Insufficient donor replication for planned exploratory contrast")
    return matched, strata


def aggregate_selected_counts(path, selected):
    """Read only selected CSR rows (hundreds), not the 610k-cell expression matrix."""
    keys = sorted(set(zip(selected.donor, selected.group, strict=True)))
    lookup = {key: i for i, key in enumerate(keys)}
    with h5py.File(path, "r") as f:
        group = f["raw/X"]
        require(group.attrs["encoding-type"] == "csr_matrix", "Raw matrix must be CSR")
        genes = read_column(f["raw/var"], f["raw/var"].attrs["_index"])
        names = read_column(f["raw/var"], "feature_name")
        current_genes = read_column(f["var"], f["var"].attrs["_index"])
        require(np.array_equal(genes, current_genes), "Raw and current gene orders differ")
        filtered = f["var/feature_is_filtered"][:]
        indptr = group["indptr"][:]
        summed = np.zeros((len(keys), len(genes)), dtype=np.int64)
        depth = []
        for cell in selected.sort_values("source_row").itertuples():
            start, stop = indptr[cell.source_row : cell.source_row + 2]
            values = group["data"][start:stop]
            indices = group["indices"][start:stop]
            require(
                np.isfinite(values).all()
                and (values >= 0).all()
                and np.array_equal(values, np.rint(values)),
                "Raw entries are not integer counts",
            )
            np.add.at(summed[lookup[(cell.donor, cell.group)]], indices, values.astype(np.int64))
            depth.append({"cell_id": cell.cell_id, "raw_total_counts": int(values.sum())})
    metadata = pd.DataFrame(keys, columns=["donor", "group"])
    metadata["pseudobulk"] = metadata.donor + "__" + metadata.group
    metadata["cells"] = [
        int(((selected.donor == d) & (selected.group == g)).sum()) for d, g in keys
    ]
    metadata["total_counts"] = summed.sum(axis=1)
    require((metadata.total_counts > 0).all(), "Empty pseudobulk")
    return (
        summed,
        metadata,
        pd.DataFrame({"gene_id": genes, "gene": names, "source_filtered": filtered}),
        pd.DataFrame(depth),
    )


def paired_de(counts, samples, genes):
    """One declared contrast; zero-variance / low-count genes are not tested."""
    cpm = counts / counts.sum(axis=1, keepdims=True) * 1e6
    logcpm = np.log2(cpm + 1)
    donors = sorted(samples.donor.unique())
    lookup = {(r.donor, r.group): r.Index for r in samples.itertuples()}
    case = np.array([logcpm[lookup[(d, "RFU_positive")]] for d in donors])
    background = np.array([logcpm[lookup[(d, "background")]] for d in donors])
    differences = case - background
    out = genes.copy()
    out["donors"] = len(donors)
    out["mean_paired_log2cpm_difference"] = differences.mean(axis=0)
    out["cpm_ge1_pseudobulks"] = (cpm >= 1).sum(axis=0)
    tcr = out.gene.str.match(r"^TR[ABDG][VDJC]")
    variance = differences.var(axis=0, ddof=1)
    out["test_status"] = np.select(
        [out.source_filtered, tcr, out.cpm_ge1_pseudobulks.lt(len(donors)), variance <= 1e-14],
        ["source_filtered", "TCR_gene_excluded", "low_expression", "zero_paired_variance"],
        default="tested",
    )
    keep = out.test_status.eq("tested").to_numpy()
    out["p_value"] = np.nan
    out["q_value"] = np.nan
    result = stats.ttest_rel(case[:, keep], background[:, keep], axis=0, alternative="two-sided")
    require(np.isfinite(result.pvalue).all(), "Nonfinite DE p-value")
    out.loc[keep, "p_value"] = result.pvalue
    out.loc[keep, "q_value"] = stats.false_discovery_control(result.pvalue, method="bh")
    out["minus_log10_q"] = -np.log10(out.q_value.clip(lower=np.finfo(float).tiny))
    out["q_lt_0_05"] = out.q_value.lt(0.05)
    out["label_gene"] = False
    ranked = out.loc[keep].assign(abs_effect=out.loc[keep, "mean_paired_log2cpm_difference"].abs())
    labels = (
        ranked.sort_values(
            ["q_value", "p_value", "abs_effect", "gene_id"], ascending=[True, True, False, True]
        )
        .head(4)
        .index
    )
    out.loc[labels, "label_gene"] = True
    return out, logcpm


def prepare(spec, out):
    verify_sources(spec)
    require(not out.exists() or not any(out.iterdir()), "Preparation directory must be empty")
    (out / "tables").mkdir(parents=True)
    (out / "manifests").mkdir()
    print("Reading frozen annotation, UMAP and assignments", flush=True)
    cells = load_cells(spec)
    qualified = cells[cells.assignment_status.eq("qualified")].copy()
    support = read_table(spec["sources"]["support"]["path"])
    ranking = rank_rfus(qualified)
    joined = ranking.merge(
        support,
        left_on="rfu",
        right_on="rfu_label",
        validate="one_to_one",
        suffixes=("", "_frozen"),
    )
    require(
        len(joined) == len(ranking)
        and joined.cells.eq(joined.cells_frozen).all()
        and joined.donors.eq(joined.donors_frozen).all(),
        "Support disagrees with frozen report",
    )
    require(
        ranking.head(3).rfu.tolist() == spec["selected_rfus"], "Metadata-only RFU selection changed"
    )
    require(
        len(cells) == 610429 and len(qualified) == 233913, "Frozen atlas/qualified totals changed"
    )
    save_table(ranking, out, "rfu_selection_ranking")
    seq, frequency, summary = logo_tables(qualified, spec["selected_rfus"])
    save_table(seq, out, "B_sequences")
    save_table(frequency, out, "B_logo_frequencies")
    save_table(summary, out, "B_rfu_support")
    seed = spec["embedding"]["seed"]
    background = set(
        stable_order(cells.cell_id.tolist(), seed)[: spec["embedding"]["background_cells"]]
    )
    cells["context_background"] = cells.cell_id.isin(background)
    cells["selected_rfu"] = cells.rfu.isin(spec["selected_rfus"])
    displayed = cells[cells.context_background | cells.selected_rfu]
    save_table(
        displayed[
            [
                "cell_id",
                "source_row",
                "umap_1",
                "umap_2",
                "cell_type",
                "donor",
                "tissue",
                "assignment_status",
                "rfu",
                "context_background",
                "selected_rfu",
            ]
        ],
        out,
        "A_umap_cells",
    )
    save_table(
        cells.groupby(["cell_type", "assignment_status"], observed=True)
        .size()
        .rename("atlas_cells")
        .reset_index(),
        out,
        "A_annotation_assignment_counts",
    )
    save_table(
        qualified[qualified.rfu.isin(spec["selected_rfus"])][
            ["cell_id", "rfu", "cdr3aa", "v_call", "donor", "tissue", "cell_type", "rfu_score"]
        ],
        out,
        "B_assigned_cells",
    )
    matched, strata = match_cells(qualified, spec["contrast"], seed)
    print(
        f"Reading raw expression for {matched.included_in_de.sum()} matched cells only", flush=True
    )
    retained = matched[matched.included_in_de].copy()
    counts, samples, genes, depth = aggregate_selected_counts(
        spec["sources"]["atlas"]["path"], retained
    )
    save_table(
        matched[
            [
                "cell_id",
                "source_row",
                "donor_id",
                "donor",
                "group",
                "rfu",
                "tissue",
                "cell_type",
                "library_id",
                "v_call",
                "stratum",
                "included_in_de",
            ]
        ].merge(depth, how="left", on="cell_id", validate="one_to_one"),
        out,
        "C_matched_cells",
    )
    save_table(strata, out, "C_matching_strata")
    save_table(samples, out, "C_pseudobulk_samples")
    de, logcpm = paired_de(counts, samples, genes)
    de["contrast"] = "RFU3526_minus_matched_other_qualified_RFUs"
    save_table(de, out, "C_differential_expression")
    for values, name in [(counts, "C_pseudobulk_counts"), (logcpm, "C_pseudobulk_log2cpm")]:
        table = pd.DataFrame(values.T, columns=samples.pseudobulk)
        table.insert(0, "gene_id", genes.gene_id)
        save_table(table, out, name)
    included_counts = samples[samples.group.eq("RFU_positive")].cells
    metrics = {
        "atlas_cells": len(cells),
        "qualified_cells": len(qualified),
        "displayed_umap_cells": len(displayed),
        "all_selected_rfu_cells": int(cells.selected_rfu.sum()),
        "de_donors": samples.donor.nunique(),
        "de_cells_per_group": int(included_counts.sum()),
        "de_cells_per_donor_min": int(included_counts.min()),
        "de_cells_per_donor_max": int(included_counts.max()),
        "de_genes_tested": int(de.test_status.eq("tested").sum()),
        "de_genes_q_lt_0_05": int(de.q_lt_0_05.sum()),
        "selected_rfus": spec["selected_rfus"],
    }
    json_save(out / "manifests/summary.json", metrics)
    json_save(
        out / "manifests/preparation.json",
        {
            "spec": spec,
            "input_hashes_verified": True,
            "code_sha256": sha256(Path(__file__)),
            "command": shlex.join(
                [
                    sys.executable,
                    "-m",
                    "manuscript.scripts.singlecell_context_figures",
                    "--prepare-only",
                    "--out",
                    str(out),
                ]
            ),
            "versions": versions(),
            "outputs": {
                str(p.relative_to(out)): sha256(p) for p in sorted(out.rglob("*")) if p.is_file()
            },
        },
    )
    print(json.dumps(metrics), flush=True)


def versions():
    return {
        name: importlib.metadata.version(name)
        for name in ["numpy", "pandas", "scipy", "h5py", "matplotlib", "Pillow", "scrfu"]
    }


def verify_prepared(root, spec):
    record = json.loads((root / "manifests/preparation.json").read_text())
    require(record["spec"] == spec, "Prepared data specification changed")
    for name, checksum in record["outputs"].items():
        require(sha256(safe_child(root, name)) == checksum, f"Changed prepared table: {name}")
    return record


def validate_source_tables(out, spec):
    """Rejoin cells to the frozen object/assignments, without expression fitting."""
    cells = load_cells(spec).set_index("cell_id")
    whole_atlas = (
        cells.groupby(["cell_type", "assignment_status"], observed=True)
        .size()
        .rename("atlas_cells")
        .reset_index()
    )
    pd.testing.assert_frame_equal(
        whole_atlas,
        read_table(out / "tables/A_annotation_assignment_counts.tsv"),
        check_dtype=False,
        check_categorical=False,
    )
    plot = read_table(out / "tables/A_umap_cells.tsv").set_index("cell_id")
    source = cells.loc[plot.index]
    for col in [
        "source_row",
        "umap_1",
        "umap_2",
        "cell_type",
        "donor",
        "tissue",
        "assignment_status",
        "rfu",
    ]:
        require(
            np.array_equal(plot[col].to_numpy(), source[col].to_numpy()),
            f"UMAP source mismatch: {col}",
        )
    require(
        set(cells.index[cells.rfu.isin(spec["selected_rfus"])]).issubset(plot.index),
        "Selected RFU cells omitted",
    )
    assigned = read_table(out / "tables/B_assigned_cells.tsv").set_index("cell_id")
    for col in assigned.columns:
        require(
            np.array_equal(assigned[col], cells.loc[assigned.index, col]),
            f"Logo cell source mismatch: {col}",
        )
    qualified = cells[cells.assignment_status.eq("qualified")].reset_index()
    for actual, name in zip(
        logo_tables(qualified, spec["selected_rfus"]),
        ["B_sequences", "B_logo_frequencies", "B_rfu_support"],
        strict=True,
    ):
        saved = read_table(out / f"tables/{name}.tsv")
        pd.testing.assert_frame_equal(
            actual, saved, check_dtype=False, check_exact=False, rtol=1e-14, atol=0
        )
    matched = read_table(out / "tables/C_matched_cells.tsv")
    for col in ["source_row", "donor", "rfu", "tissue", "cell_type", "library_id", "v_call"]:
        require(
            np.array_equal(matched[col], cells.loc[matched.cell_id, col]),
            f"Matched cell source mismatch: {col}",
        )
    retained = matched[matched.included_in_de]
    require(retained.cell_id.is_unique, "Repeated matched cell")
    sizes = retained.groupby(["stratum", "group"], observed=True).size().unstack(fill_value=0)
    require(sizes.RFU_positive.eq(sizes.background).all(), "Unbalanced matching stratum")
    require(
        retained[retained.group.eq("RFU_positive")].rfu.eq(spec["contrast"]["case_rfu"]).all(),
        "Case RFU mismatch",
    )
    require(
        retained[retained.group.eq("background")].rfu.ne(spec["contrast"]["case_rfu"]).all(),
        "RFU-positive background",
    )
    counts, samples, genes, _ = aggregate_selected_counts(
        spec["sources"]["atlas"]["path"], retained
    )
    recorded = read_table(out / "tables/C_pseudobulk_counts.tsv")
    require(
        np.array_equal(recorded.gene_id, genes.gene_id)
        and np.array_equal(recorded.iloc[:, 1:].to_numpy().T, counts),
        "Pseudobulk counts differ from raw source rows",
    )
    pd.testing.assert_frame_equal(
        samples, read_table(out / "tables/C_pseudobulk_samples.tsv"), check_dtype=False
    )
    # Verify transformation, effect sizes and BH arithmetic; do not repeat tests/models.
    logcpm = np.log2(counts / counts.sum(axis=1, keepdims=True) * 1e6 + 1)
    saved_log = read_table(out / "tables/C_pseudobulk_log2cpm.tsv").iloc[:, 1:].to_numpy().T
    require(np.allclose(logcpm, saved_log, atol=0, rtol=1e-14), "Log-CPM transform mismatch")
    lookup = samples.set_index(["donor", "group"]).pseudobulk
    effect = np.mean(
        [
            saved_log[samples.pseudobulk.eq(lookup[d, "RFU_positive"])][0]
            - saved_log[samples.pseudobulk.eq(lookup[d, "background"])][0]
            for d in sorted(samples.donor.unique())
        ],
        axis=0,
    )
    de = pd.read_csv(
        out / "tables/C_differential_expression.tsv", sep="\t", float_precision="round_trip"
    )
    require(
        np.allclose(de.mean_paired_log2cpm_difference, effect, atol=1e-14),
        "DE effects differ from frozen pseudobulks",
    )
    tested = de[de.test_status.eq("tested")]
    summary = json.loads((out / "manifests/summary.json").read_text())
    donor_counts = samples[samples.group.eq("RFU_positive")].cells
    expected_summary = {
        "atlas_cells": len(cells),
        "qualified_cells": int(cells.assignment_status.eq("qualified").sum()),
        "displayed_umap_cells": len(plot),
        "all_selected_rfu_cells": len(assigned),
        "de_donors": samples.donor.nunique(),
        "de_cells_per_group": int(donor_counts.sum()),
        "de_cells_per_donor_min": int(donor_counts.min()),
        "de_cells_per_donor_max": int(donor_counts.max()),
        "de_genes_tested": len(tested),
        "de_genes_q_lt_0_05": int(tested.q_value.lt(0.05).sum()),
        "selected_rfus": spec["selected_rfus"],
    }
    require(summary == expected_summary, "Displayed summary differs from source tables")
    # Independent Benjamini-Hochberg check on the stored p-values.
    order = np.argsort(tested.p_value.to_numpy(), kind="stable")
    expected = np.minimum.accumulate(
        (tested.p_value.to_numpy()[order] * len(tested) / np.arange(1, len(tested) + 1))[::-1]
    )[::-1].clip(max=1)
    require(
        np.allclose(tested.q_value.to_numpy()[order], expected, atol=1e-15),
        "BH adjustment mismatch",
    )
    return {
        "frozen_umap_coordinates_exact": len(plot),
        "selected_rfu_cells_verified": len(assigned),
        "raw_expression_rows_verified": len(retained),
        "tested_genes_verified": len(tested),
        "assignment_or_embedding_recomputed": False,
        "status": "passed",
    }


def render(spec, prepared, out):
    from ._singlecell_context_drawing import render_all

    verify_sources(spec)
    verify_prepared(prepared, spec)
    require(
        not out.exists() or not any(out.iterdir()),
        "Render output must be empty; use --verify-only for completed packs",
    )
    out.mkdir(parents=True, exist_ok=True)
    shutil.copytree(prepared / "tables", out / "tables")
    shutil.copytree(prepared / "manifests", out / "manifests")
    summary = json.loads((out / "manifests/summary.json").read_text())
    validation = validate_source_tables(out, spec)
    exports = render_all(out, spec, summary)
    shutil.copyfile(DOC, out / "manifests/captions.md")
    shutil.copyfile(REPO / "manuscript/figures/singlecell_context_README.md", out / "README.md")
    (out / "tables/README.md").write_text(
        "# Source tables\n\nSee ../README.md for table units, denominators and column descriptions, "
        "../manifests/table_columns.json for exact headers, and ../manifests/preparation.json "
        "for input hashes and the preparation command. Tables are external biological data; "
        "do not commit them to Git.\n"
    )
    schema = {p.name: list(read_table(p).columns) for p in sorted((out / "tables").glob("*.tsv"))}
    json_save(out / "manifests/table_columns.json", schema)
    command = shlex.join(
        [
            sys.executable,
            "-m",
            "manuscript.scripts.singlecell_context_figures",
            "--prepared",
            str(prepared),
            "--out",
            str(out),
        ]
    )
    json_save(
        out / "manifests/figure_manifest.json",
        {
            "spec": spec,
            "exports": exports,
            "summary": summary,
            "command": command,
            "preparation_manifest": "manifests/preparation.json",
            "versions": versions(),
            "source_validation": validation,
            "code_sha256": {
                str(p.relative_to(REPO)): sha256(p)
                for p in [
                    Path(__file__),
                    Path(__file__).with_name("_singlecell_context_drawing.py"),
                    Path(__file__).with_name("_rp1_14_figure_drawing.py"),
                    Path(__file__).with_name("rp1_14_reusable_figures.py"),
                    SPEC,
                    DOC,
                    REPO / "manuscript/figures/singlecell_context_README.md",
                ]
            },
        },
    )
    json_save(
        out / "completion.json",
        {
            "status": "complete",
            "outputs": {
                str(p.relative_to(out)): sha256(p) for p in sorted(out.rglob("*")) if p.is_file()
            },
        },
    )
    print(json.dumps(verify_pack(out, spec)), flush=True)


def verify_pack(out, spec):
    from ._singlecell_context_drawing import verify_exports

    completed = json.loads((out / "completion.json").read_text())
    for name, checksum in completed["outputs"].items():
        require(sha256(safe_child(out, name)) == checksum, f"Output hash mismatch: {name}")
    require(completed["status"] == "complete", "Incomplete pack")
    verify_prepared(out, spec)
    manifest = json.loads((out / "manifests/figure_manifest.json").read_text())
    verify_exports(out, manifest["exports"])
    require(manifest["spec"] == spec, "Rendering specification drift")
    return {
        "status": "passed",
        "exports": len(manifest["exports"]),
        "outputs_hashed": len(completed["outputs"]),
        "source_validation": manifest["source_validation"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path)
    parser.add_argument(
        "--prepared", type=Path, help="Frozen figure-specific tables; rendering does no DE"
    )
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    spec = json.loads(SPEC.read_text())
    out = (args.out or Path(spec["output"])).resolve()
    require(not out.is_relative_to(REPO), "Biological tables/figures must remain outside Git")
    require(not (args.prepare_only and (args.prepared or args.verify_only)), "Conflicting modes")
    if args.verify_only:
        verify_sources(spec)
        print(json.dumps(verify_pack(out, spec)))
    elif args.prepare_only:
        prepare(spec, out)
    else:
        require(
            args.prepared is not None, "First freeze the bounded preparation with --prepare-only"
        )
        render(spec, args.prepared.resolve(), out)


if __name__ == "__main__":
    main()
