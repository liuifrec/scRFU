"""Frozen external repertoire application; no cell-state inference or tuning.

Run assignment and measurement as separate restartable stages. All biological
outputs stay external. Required-visit support is always evaluated on threshold
cells, including for nearest-label sensitivity. Missing is never encoded as zero.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import inspect
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import scrfu
from scrfu.io import file_sha256
from scrfu.tl import (
    assign_rfu,
    multiscale_repertoire_change,
    permute_fixed_groups,
    validate_completed_rfu_run,
)

from ._figure_common import clean_axis, panel_label, save_figure, set_style
from .radiation_methods import REPO, assert_hash, completed, json_save, save, validate_boundary

FROZEN_CONFIG = "cb51b5933bb1494401934f016f8f964ef4327e62f842d3d931d575587fa89e8b"
INTERVALS = {
    (1, 2): "pre_to_last_radiation",
    (2, 3): "last_radiation_to_6weeks",
    (1, 3): "pre_to_6weeks",
}


def frozen_config(prepared: Path) -> dict:
    path = prepared / "frozen_development_configuration.json"
    assert_hash(path, FROZEN_CONFIG)
    config = json.loads(path.read_text())
    if config != json.loads((REPO / "manuscript/config/radiation_methods_v1.json").read_text()):
        raise ValueError("External settings differ from the frozen development configuration.")
    return config


def verified_stage(path: Path) -> dict:
    manifest = json.loads((path / "completion.json").read_text())
    if not completed(path, manifest["fingerprint"]):
        raise ValueError("Input stage is not complete.")
    return manifest


def environment() -> dict:
    return {
        "python": sys.executable,
        "scrfu_import": scrfu.__file__,
        "git_sha": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
        ).strip(),
        "R": subprocess.check_output(
            ["/usr/bin/Rscript", "--version"], stderr=subprocess.STDOUT, text=True
        ).strip(),
        "packages": {
            x: importlib.metadata.version(x)
            for x in ("numpy", "pandas", "scrfu", "matplotlib", "pyarrow")
        },
    }


def identity(inputs: dict[str, Path], functions: tuple) -> dict:
    # Assignment fingerprint excludes downstream plotting/measurement edits.
    return {
        "inputs": {
            k: {"path": str(v.resolve()), "sha256": file_sha256(v)} for k, v in inputs.items()
        },
        "functions": {
            f.__name__: hashlib.sha256(inspect.getsource(f).encode()).hexdigest() for f in functions
        },
        "core_code": {
            str(p.relative_to(REPO)): file_sha256(p)
            for p in sorted((REPO / "src/scrfu").rglob("*.py"))
        },
        "packages": environment()["packages"],
    }


def start(out: Path, scientific: dict) -> str | None:
    if out.resolve().is_relative_to(REPO):
        raise ValueError("Biological results must remain outside Git.")
    fingerprint = hashlib.sha256(json.dumps(scientific, sort_keys=True).encode()).hexdigest()
    if completed(out, fingerprint):
        print(f"Verified completed {out.name}; no recomputation.", flush=True)
        return None
    out.mkdir(parents=True, exist_ok=True)
    json_save({"status": "running", "fingerprint": fingerprint}, out / "completion.json")
    json_save(
        {
            "scientific_identity": scientific,
            "environment": environment(),
            "command": sys.argv,
            "script_sha256": file_sha256(Path(__file__)),
        },
        out / "provenance.json",
    )
    return fingerprint


def finish(out: Path, fingerprint: str) -> None:
    # Backend trees are validated by the core validator; top-level deliverables
    # have independent hashes and are sealed only after all exports finish.
    outputs = {
        p.name: file_sha256(p)
        for p in sorted(out.iterdir())
        if p.is_file() and p.name != "completion.json"
    }
    json_save(
        {"status": "complete", "fingerprint": fingerprint, "outputs": outputs},
        out / "completion.json",
    )


def check_cells(frame: pd.DataFrame, registry: pd.DataFrame) -> None:
    validate_boundary(frame, registry)
    if set(frame.dataset_id) != {"GSE280982"} or not set(frame.compartment) <= {"tumor", "blood"}:
        raise ValueError("External dataset/compartment boundary violated.")
    if frame.cell_id.duplicated().any() or frame.input_row_id.duplicated().any():
        raise ValueError("Primary receptor/cell IDs must be unique.")
    if not frame.cell_id.eq(frame.sample_id + ":" + frame.source_barcode).all():
        raise ValueError("Sample-barcode namespace mismatch.")
    lookup = registry.set_index("sample_id")
    if lookup.index.duplicated().any():
        raise ValueError("Ambiguous sample registry.")
    for key in ("donor", "visit", "compartment"):
        if not frame[key].eq(frame.sample_id.map(lookup[key])).all():
            raise ValueError("Sample/visit metadata disagree.")


def assign(prepared: Path, out: Path, reference: Path) -> None:
    verified_stage(prepared)
    config = frozen_config(prepared)
    sources = {
        "prepared": prepared / "completion.json",
        "config": prepared / "frozen_development_configuration.json",
    }
    for name, key in (
        ("RFU.R", "rfu_r_sha256"),
        ("km5000noMax.Rdata", "km5000_rdata_sha256"),
        ("trimerMDSfit_small.Rdata", "trimer_rdata_sha256"),
    ):
        sources[name] = reference / name
        assert_hash(sources[name], config["reference"][key])
    scientific = identity(sources, (assign, frozen_config, check_cells))
    fingerprint = start(out, scientific)
    if fingerprint is None:
        validate_completed_rfu_run(out, expected_provenance=config["reference"])
        return
    frame = pd.read_parquet(prepared / "gse280982_primary_trb_matched_gex.parquet")
    registry = pd.read_csv(prepared / "gse280982_verified_sample_registry.tsv", sep="\t")
    check_cells(frame, registry)
    result = assign_rfu(
        frame,
        rfu_dir=reference,
        mode="standard",
        threshold=0.6,
        deduplicate=True,
        chunk_size=500,
        max_workers=1,
        resume=True,
        workdir=out / "backend",
        wrapper_r_path=REPO / "r/run_rfu_repo.R",
        rscript_bin="/usr/bin/Rscript",
    )
    for name, table in (
        ("receptors", frame),
        ("unique_sequence_map", result.mapping),
        ("rfu_results_per_sequence", result.per_sequence),
        ("rfu_results_per_row", result.per_row),
    ):
        table.to_csv(out / f"{name}.tsv.gz", sep="\t", index=False)
    json_save(
        {**result.provenance, "workflow": "frozen_GSE280982_external"}, out / "run_manifest.json"
    )
    json_save(
        validate_completed_rfu_run(out, expected_provenance=config["reference"]),
        out / "validation.json",
    )
    finish(out, fingerprint)
    print("External assignments complete and independently validated.", flush=True)


def attach(frame: pd.DataFrame, rows: pd.DataFrame) -> pd.DataFrame:
    fields = [
        "input_row_id",
        "rfu_label_nearest",
        "rfu_score",
        "rfu_pass_threshold",
        "eligibility_status",
    ]
    result = frame.merge(rows[fields], on="input_row_id", validate="one_to_one", how="left")
    if result.eligibility_status.isna().any() or len(result) != len(rows):
        raise ValueError("Assignment/receptor coverage mismatch.")
    result["rfu_pass_threshold"] = result.rfu_pass_threshold.astype("boolean").fillna(False)
    if result[["junction", "v_call", "j_call"]].isna().any().any():
        raise ValueError("Source nucleotide/V/J identity incomplete; do not guess clones.")
    result["clone"] = result.junction + "|" + result.v_call + "|" + result.j_call
    result["rfu"] = result.rfu_label_nearest
    result["vj"] = result.v_call + "|" + result.j_call
    result["length_bin"] = result.cdr3aa.str.len() // 5
    known = result[result.rfu.notna()]
    if known.groupby("clone")[["rfu", "rfu_pass_threshold"]].nunique().gt(1).any().any():
        raise ValueError("Receptor group or threshold eligibility varies across visits.")
    return result


def support_tables(
    frame: pd.DataFrame, registry: pd.DataFrame, minimum: int = 100
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if minimum != 100:
        raise ValueError("Do not lower or tune the frozen 100-qualified-cell support rule.")
    rows = []
    for (donor, compartment), part in registry.groupby(["donor", "compartment"], sort=True):
        for visit in (1, 2, 3):
            source = part[part.visit.eq(visit)]
            observed = frame[
                frame.donor.eq(donor) & frame.compartment.eq(compartment) & frame.visit.eq(visit)
            ]
            available = bool(
                len(source) and source.status.iloc[0] == "prepared_before_RFU_assignment"
            )
            if available and len(observed) != source.primary_TRB_cells_in_GEX.iloc[0]:
                raise ValueError("Prepared sample cell denominator changed.")
            qualified = observed[observed.rfu_pass_threshold]
            rows.append(
                {
                    "donor": donor,
                    "compartment": compartment,
                    "visit": visit,
                    "source_available": available,
                    "primary_TRB_cells": len(observed) if available else np.nan,
                    "qualified_cells": len(qualified) if available else np.nan,
                    "coverage": len(qualified) / len(observed)
                    if available and len(observed)
                    else np.nan,
                    "gex_cells": source.gex_barcodes.iloc[0] if available else np.nan,
                    "qualified_receptors": qualified.clone.nunique() if available else np.nan,
                    "qualified_RFUs": qualified.rfu.nunique() if available else np.nan,
                    "status": "supported"
                    if available and len(qualified) >= minimum
                    else "below_frozen_support"
                    if available
                    else "unavailable",
                }
            )
    coverage = pd.DataFrame(rows)
    pairs = []
    for (donor, compartment), group in coverage.groupby(["donor", "compartment"], sort=True):
        group = group.set_index("visit")
        for (a, b), interval in INTERVALS.items():
            statuses = group.loc[[a, b], "status"]
            pairs.append(
                {
                    "donor": donor,
                    "compartment": compartment,
                    "visit_before": a,
                    "visit_after": b,
                    "interval": interval,
                    "before_qualified_cells": group.loc[a, "qualified_cells"],
                    "after_qualified_cells": group.loc[b, "qualified_cells"],
                    "status": "analyzed"
                    if statuses.eq("supported").all()
                    else "unavailable_required_visit"
                    if statuses.eq("unavailable").any()
                    else "insufficient_qualified_cells",
                }
            )
    return coverage, pd.DataFrame(pairs)


def counts(frame: pd.DataFrame, visit: int, weighting: str) -> pd.Series:
    x = frame[frame.visit.eq(visit)].groupby("clone").size()
    return x if weighting == "cell" else x.clip(upper=1)


def persistence(person: pd.DataFrame, pair: dict, detections: list[int]) -> list[dict]:
    result = []
    for rfu, group in person.groupby("rfu", sort=True):
        a, b = (counts(group, pair[v], "cell") for v in ("visit_before", "visit_after"))
        for detection in detections:
            result.append(
                {
                    **pair,
                    "rfu": rfu,
                    "detection_min_cells": detection,
                    "before_cells": int(a.sum()),
                    "after_cells": int(b.sum()),
                    "observed_before": bool(a.sum() >= detection),
                    "observed_after": bool(b.sum() >= detection),
                    "observed_shared_receptors": len(a.index.intersection(b.index)),
                    "dominant_fraction_before": a.max() / a.sum() if len(a) else np.nan,
                    "dominant_fraction_after": b.max() / b.sum() if len(b) else np.nan,
                }
            )
    return result


def measurements(frame: pd.DataFrame, pairs: pd.DataFrame, config: dict, out: Path) -> None:
    eligible_pairs = pairs[pairs.status.eq("analyzed")].drop(columns="status")
    mapping = frame[frame.rfu.notna()].drop_duplicates("clone").set_index("clone")
    primary = frame[frame.rfu_pass_threshold]
    qualified_map = primary.drop_duplicates("clone").set_index("clone")
    strata = qualified_map[["v_call", "j_call", "length_bin"]]
    variation = (
        strata.assign(rfu=qualified_map.rfu).groupby(list(strata.columns)).rfu.transform("nunique")
    )
    random_maps = [
        permute_fixed_groups(
            qualified_map.rfu, strata=strata, random_state=config["random_seed"] + i
        )
        for i in range(config["random_group_replicates"])
    ]
    distances, persistent, similarities, controls, depths, dominance = [], [], [], [], [], []
    for pair_i, pair in enumerate(eligible_pairs.to_dict("records")):
        subset = frame[
            frame.donor.eq(pair["donor"])
            & frame.compartment.eq(pair["compartment"])
            & frame.visit.isin([pair["visit_before"], pair["visit_after"]])
        ]
        for policy in ("threshold", "nearest"):
            person = (
                subset[subset.rfu_pass_threshold]
                if policy == "threshold"
                else subset[subset.rfu.notna()]
            )
            for weight in config["weights"]:
                a, b = (counts(person, pair[v], weight) for v in ("visit_before", "visit_after"))
                base = {**pair, "policy": policy, "weighting": weight}
                for column, label in (("rfu", "RFU"), ("v_call", "TRBV"), ("vj", "TRBV_TRBJ")):
                    distances.append(
                        {
                            **base,
                            "grouping": label,
                            **multiscale_repertoire_change(a, b, mapping[column]).summary,
                        }
                    )
                aa, bb = (
                    a.groupby(mapping.rfu).sum().align(b.groupby(mapping.rfu).sum(), fill_value=0)
                )
                similarities.append(
                    {
                        **base,
                        "rfu_cosine": float(
                            np.dot(aa, bb) / (np.linalg.norm(aa) * np.linalg.norm(bb))
                        ),
                    }
                )
            if policy == "threshold":
                persistent.extend(persistence(person, pair, config["rfu_detection_cell_counts"]))
        person = subset[subset.rfu_pass_threshold]
        a, b = (counts(person, pair[v], "cell") for v in ("visit_before", "visit_after"))
        for i, group_map in enumerate(random_maps):
            observed = a.index.union(b.index)
            controls.append(
                {
                    **pair,
                    "replicate": i,
                    "seed": config["random_seed"] + i,
                    "global_map_exchangeable_fraction": float(variation.gt(1).mean()),
                    "pair_receptors_exchangeable_fraction": float(
                        variation.loc[observed].gt(1).mean()
                    ),
                    "pair_labels_changed_fraction": float(
                        group_map.loc[observed].ne(qualified_map.loc[observed, "rfu"]).mean()
                    ),
                    **multiscale_repertoire_change(a, b, group_map).summary,
                }
            )
        cells_a, cells_b = (
            person.loc[person.visit.eq(pair[v]), "clone"].to_numpy()
            for v in ("visit_before", "visit_after")
        )
        n = min(config["depth_cap_cells"], len(cells_a), len(cells_b))
        for i in range(config["depth_replicates"]):
            seed = config["random_seed"] + 10000 * (pair_i + 1) + i
            rng = np.random.default_rng(seed)
            x, y = (
                pd.Series(rng.choice(c, n, replace=False)).value_counts()
                for c in (cells_a, cells_b)
            )
            depths.append(
                {
                    **pair,
                    "replicate": i,
                    "seed": seed,
                    "cells_per_visit": n,
                    "scheme": "observed_cell_subsampling_without_replacement",
                    **multiscale_repertoire_change(x, y, mapping.rfu).summary,
                }
            )
        remove = {a.idxmax(), b.idxmax()}
        x, y = (s.drop(list(remove), errors="ignore") for s in (a, b))
        dominance.append(
            {
                **pair,
                "scheme": "remove_union_of_each_visit_dominant_primary_TRB",
                "removed_receptors": len(remove),
                "dominant_fraction_before": float(a.max() / a.sum()),
                "dominant_fraction_after": float(b.max() / b.sum()),
                "retained_fraction_before": float(x.sum() / a.sum()),
                "retained_fraction_after": float(y.sum() / b.sum()),
                **multiscale_repertoire_change(x, y, mapping.rfu).summary,
            }
        )
    for name, rows in (
        ("multiscale_pairs", distances),
        ("rfu_persistence", persistent),
        ("similarity", similarities),
        ("fixed_group_controls", controls),
        ("cell_depth_sensitivity", depths),
        ("clone_dominance_sensitivity", dominance),
    ):
        save(pd.DataFrame(rows), out, name + ".tsv")
    summary = []
    p = pd.DataFrame(persistent)
    keys = ["donor", "compartment", "interval", "detection_min_cells"]
    if len(p):
        for key, group in p.groupby(keys):
            union = group.observed_before | group.observed_after
            shared = group.observed_before & group.observed_after
            summary.append(
                {
                    **dict(zip(keys, key, strict=True)),
                    "detected_union_RFUs": int(union.sum()),
                    "persistent_RFUs": int(shared.sum()),
                    "persistent_over_union": shared.sum() / union.sum() if union.any() else np.nan,
                    "persistent_without_shared_receptor": int(
                        (shared & group.observed_shared_receptors.eq(0)).sum()
                    ),
                    "persistent_without_shared_receptor_fraction": group.loc[
                        shared, "observed_shared_receptors"
                    ]
                    .eq(0)
                    .mean(),
                }
            )
    save(pd.DataFrame(summary), out, "persistence_summary.tsv")


def external_figure(out: Path) -> None:
    import matplotlib.pyplot as plt

    set_style()
    coverage = pd.read_csv(out / "visit_support.tsv", sep="\t")
    d = pd.read_csv(out / "multiscale_pairs.tsv", sep="\t")
    d = d[d.policy.eq("threshold") & d.weighting.eq("cell") & d.grouping.eq("RFU")].reset_index(
        drop=True
    )
    p = pd.read_csv(out / "persistence_summary.tsv", sep="\t")
    p = p[p.detection_min_cells.eq(1)]
    depth = pd.read_csv(out / "cell_depth_sensitivity.tsv", sep="\t")
    control = pd.read_csv(out / "fixed_group_controls.tsv", sep="\t")
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    ax = axes[0, 0]
    for (donor, tissue), rows in coverage.groupby(["donor", "compartment"]):
        ax.plot(
            rows.visit,
            rows.qualified_cells,
            marker="o",
            linestyle="-" if tissue == "tumor" else "--",
            label=f"{donor} {tissue}",
        )
    ax.axhline(100, color="black", linestyle=":", linewidth=1)
    ax.set(
        xticks=[1, 2, 3],
        xticklabels=["Pre", "Last RT day", "~6 weeks"],
        ylabel="Qualified primary-TRB cells",
        title="Visit eligibility (100 cells)",
    )
    ax.legend(ncol=2, fontsize=6)
    ax = axes[0, 1]
    labels = [
        f"{r.donor} {r.compartment}\n{int(r.visit_before)}→{int(r.visit_after)}"
        for r in d.itertuples()
    ]
    x = np.arange(len(d))
    ax.plot(x, d.d_clone, "o", label="Receptor TV")
    ax.plot(x, d.d_group, "s", label="RFU TV")
    for i, r in d.iterrows():
        ax.plot([i, i], [r.d_group, r.d_clone], color="gray", linewidth=1)
    ax.set(
        xticks=x,
        xticklabels=labels,
        ylabel="Total variation",
        ylim=(0, 1.05),
        title="Separate donor/interval measurements",
    )
    ax.tick_params(axis="x", labelrotation=50)
    ax.legend()
    ax = axes[1, 0]
    for row in d.itertuples():
        q = p[
            p.donor.eq(row.donor) & p.compartment.eq(row.compartment) & p.interval.eq(row.interval)
        ].iloc[0]
        ax.scatter(
            q.persistent_over_union,
            q.persistent_without_shared_receptor_fraction,
            marker="s" if row.compartment == "blood" else "o",
        )
        ax.annotate(
            f"{row.donor} {int(row.visit_before)}→{int(row.visit_after)}",
            (q.persistent_over_union, q.persistent_without_shared_receptor_fraction),
            fontsize=6,
            xytext=(3, 3),
            textcoords="offset points",
        )
    ax.set(
        xlabel="Persistent RFUs / detected union (≥1 cell)",
        ylabel="Fraction of persistent RFUs without\na shared observed receptor",
        xlim=(0, 1),
        ylim=(0, 1.05),
        title="Persistence is not receptor preservation",
    )
    ax = axes[1, 1]
    for i, r in d.iterrows():
        dep = depth[
            depth.donor.eq(r.donor)
            & depth.compartment.eq(r.compartment)
            & depth.interval.eq(r.interval)
        ]
        ctrl = control[
            control.donor.eq(r.donor)
            & control.compartment.eq(r.compartment)
            & control.interval.eq(r.interval)
        ]
        ax.plot(
            i - 0.12,
            r.aggregation_cancellation,
            "o",
            color="#0072B2",
            label="Observed cancellation" if i == 0 else None,
        )
        ax.plot(
            i + 0.12,
            ctrl.aggregation_cancellation.median(),
            "s",
            color="gray",
            label="Matched-map median" if i == 0 else None,
        )
        ax.plot(
            i,
            dep.aggregation_cancellation.median(),
            "^",
            color="#D89000",
            label="Matched-depth median" if i == 0 else None,
        )
    ax.set(
        xticks=x,
        xticklabels=labels,
        ylabel="Aggregation cancellation",
        title="Controls and empirical cell subsampling",
    )
    ax.tick_params(axis="x", labelrotation=50)
    ax.legend()
    for letter, ax in zip("ABCD", axes.flat, strict=True):
        panel_label(ax, letter)
        clean_axis(ax)
    save_figure(fig, out, "gse280982_external")
    plt.close(fig)
    sources = {
        "A": "visit_support.tsv",
        "B": "multiscale_pairs.tsv",
        "C": "persistence_summary.tsv",
        "D": "fixed_group_controls.tsv;cell_depth_sensitivity.tsv;multiscale_pairs.tsv",
    }
    save(
        pd.DataFrame(
            [
                {
                    "panel": k,
                    "figure": "gse280982_external.pdf",
                    "source_tables": v,
                    "source_hashes": json.dumps(
                        {n: file_sha256(out / n) for n in v.split(";")}, sort_keys=True
                    ),
                    "unit": "donor trajectory; visits/intervals/resamples not independent participants",
                    "limit": "qualified primary TRB; relative sampled frequencies; no antigen/causal interpretation",
                    "generating_command": "python -m manuscript.scripts.gse280982_analysis measure --workspace WORKSPACE",
                    "provenance": "provenance.json; completion.json",
                }
                for k, v in sources.items()
            ]
        ),
        out,
        "figure_source_index.tsv",
    )


def measure(prepared: Path, assignment: Path, out: Path) -> None:
    verified_stage(prepared)
    verified_stage(assignment)
    config = frozen_config(prepared)
    validate_completed_rfu_run(assignment, expected_provenance=config["reference"])
    scientific = identity(
        {
            "prepared": prepared / "completion.json",
            "assignment": assignment / "completion.json",
            "config": prepared / "frozen_development_configuration.json",
        },
        (
            measure,
            attach,
            check_cells,
            support_tables,
            counts,
            persistence,
            measurements,
            external_figure,
        ),
    )
    fingerprint = start(out, scientific)
    if fingerprint is None:
        return
    frame = pd.read_parquet(prepared / "gse280982_primary_trb_matched_gex.parquet")
    registry = pd.read_csv(prepared / "gse280982_verified_sample_registry.tsv", sep="\t")
    check_cells(frame, registry)
    frame = attach(frame, pd.read_csv(assignment / "rfu_results_per_row.tsv.gz", sep="\t"))
    coverage, pairs = support_tables(frame, registry, config["min_cells_per_visit"])
    save(coverage, out, "visit_support.tsv")
    save(pairs, out, "pair_support.tsv")
    measurements(frame, pairs, config, out)
    d = pd.read_csv(out / "multiscale_pairs.tsv", sep="\t")
    if not d.d_group.le(d.d_clone + 1e-12).all():
        raise ValueError("TV contraction failed.")
    json_save(
        {
            "primary_TRB_cells": len(frame),
            "qualified_cells": int(frame.rfu_pass_threshold.sum()),
            "qualified_RFUs": frame.loc[frame.rfu_pass_threshold, "rfu"].nunique(),
            "supported_visits": int(coverage.status.eq("supported").sum()),
            "analyzed_pairs": int(pairs.status.eq("analyzed").sum()),
            "measurement_rows": len(d),
            "all_contraction_rows_pass": True,
            "author_cell_states": "not_in_prepared_sources; see metadata_availability.json",
            "interpretation": "transportability; no population effect estimate; no absolute depletion",
        },
        out / "evidence_counts.json",
    )
    external_figure(out)
    save(
        pd.DataFrame(
            [
                {
                    "file": p.name,
                    "sha256": file_sha256(p),
                    "unit": "see column names; donor is biological replicate",
                }
                for p in sorted(out.glob("*.tsv"))
            ]
        ),
        out,
        "source_table_index.tsv",
    )
    finish(out, fingerprint)
    print("External measurements, figures and source tables complete.", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["assign", "measure"])
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--reference", type=Path, default=Path("/home/liuyuchen/ext/RFU-official"))
    args = parser.parse_args()
    prepared = args.workspace / "prepared/gse280982_v1"
    assignment = args.workspace / "prepared/gse280982_rfu_v1"
    if args.stage == "assign":
        assign(prepared, assignment, args.reference)
    else:
        measure(prepared, assignment, args.workspace / "results/gse280982_v1")


if __name__ == "__main__":
    main()
