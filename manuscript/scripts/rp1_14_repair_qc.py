"""Supplementary alignment audit and stratified bounded parity, not new endpoints.

The sampling design is independent of current-reference assay results: 1,600
unique canonical productive CDR3aa, round-robin over occupied repair-status ×
sample (donor/visit/compartment) × score × RFU-frequency strata, seed 20260916.
All sampled receptors are retained in the reported denominator, including any
mismatches. Historical assignments and biological analyses are never overwritten.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scrfu.io import file_sha256
from scrfu.tl import assign_rfu, validate_completed_rfu_run

from .gse280982_analysis import frozen_config, identity, start, verified_stage
from .radiation_methods import REPO, assert_hash, json_save, save, validate_boundary

STRATA = ["donor", "visit", "compartment", "repair_status", "score_bin", "rfu_frequency_bin"]


def vector_hash(values: pd.Series) -> str:
    array = pd.to_numeric(values).to_numpy(dtype="float64", na_value=np.nan).astype("<f8")
    array[np.isnan(array)] = np.nan
    return hashlib.sha256(array.tobytes()).hexdigest()


def vector_audit(saved: pd.DataFrame, repaired: pd.DataFrame) -> dict:
    if len(saved) != len(repaired):
        raise ValueError("Export/source vector length changed.")
    result = {"rows": len(saved)}
    for original, current in [("rfu", "rfu_numeric"), ("max_cor", "max_cor")]:
        a, b = vector_hash(saved[original]), vector_hash(repaired[current])
        result[f"{original}_original_sha256"] = a
        result[f"{original}_repaired_sha256"] = b
        if a != b:
            raise ValueError("Historical RFU or score vector was altered.")
    return result


def conflict_count(frame: pd.DataFrame, sequence: str, value: str) -> int:
    return int(frame.groupby(sequence)[value].nunique().gt(1).sum())


def stratified_sample(
    frame: pd.DataFrame, target: int = 1600, seed: int = 20260916
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    eligible = frame[frame.productive_canonical & frame.rfu.notna() & frame.max_cor.notna()].copy()
    freq = eligible.drop_duplicates("aminoAcid").groupby("rfu").size()
    low, high = freq.quantile([0.25, 0.75]).tolist()
    eligible["rfu_unique_aa_frequency"] = eligible.rfu.map(freq)
    eligible["rfu_frequency_bin"] = np.where(
        eligible.rfu_unique_aa_frequency.le(low),
        "rare_Q1",
        np.where(eligible.rfu_unique_aa_frequency.gt(high), "common_Q4", "middle"),
    )
    eligible["score_bin"] = pd.cut(
        eligible.max_cor,
        [-np.inf, 0.55, 0.59, 0.60, 0.61, 0.65, np.inf],
        right=False,
        labels=[
            "below_0.55",
            "0.55_to_0.59",
            "0.59_to_0.60",
            "0.60_to_0.61",
            "0.61_to_0.65",
            "at_least_0.65",
        ],
    ).astype(str)
    eligible = eligible.sort_values(STRATA + ["aminoAcid", "source_row"]).reset_index(drop=True)
    rng = np.random.default_rng(seed)
    groups = [
        list(rng.permutation(g.index.to_numpy())) for _, g in eligible.groupby(STRATA, sort=True)
    ]
    order = rng.permutation(len(groups))
    chosen, used = [], set()
    while len(chosen) < min(target, eligible.aminoAcid.nunique()):
        added = False
        for j in order:
            while groups[j]:
                i = groups[j].pop()
                aa = eligible.at[i, "aminoAcid"]
                if aa not in used:
                    chosen.append(i)
                    used.add(aa)
                    added = True
                    break
            if len(chosen) == target:
                break
        if not added:
            break
    sample = eligible.loc[chosen].reset_index(drop=True)
    strata = (
        eligible.groupby(STRATA, observed=True)
        .aminoAcid.nunique()
        .rename("eligible_unique_aa")
        .to_frame()
    )
    strata = (
        strata.join(sample.groupby(STRATA, observed=True).size().rename("selected_unique_aa"))
        .fillna(0)
        .reset_index()
    )
    return (
        sample,
        strata,
        {
            "target": target,
            "seed": seed,
            "eligible_unique_aa": eligible.aminoAcid.nunique(),
            "selected_unique_aa": len(sample),
            "rfu_unique_aa_Q25": low,
            "rfu_unique_aa_Q75": high,
            "occupied_strata": len(strata),
            "sampled_strata": int(strata.selected_unique_aa.gt(0).sum()),
            "stratum_attribution": "one source occurrence per selected AA; not independent assays for repeated occurrences",
            "sampling": "fixed-seed round-robin across occupied strata; no current score/label used",
        },
    )


def run(workspace: Path, reference: Path, historical: Path) -> None:
    prepared = workspace / "prepared/rp1_14_v1"
    verified_stage(prepared)
    config = frozen_config(workspace / "prepared/gse280982_v1")
    original_maps = {
        p.name.replace(".cdr3_rfu.tsv", ""): p
        for p in (REPO / "data/RFU_out_RP_P1-14").rglob("*.cdr3_rfu.tsv")
    }
    inputs = {
        "prepared": prepared / "completion.json",
        "historical_code": historical,
        "config": workspace / "prepared/gse280982_v1/frozen_development_configuration.json",
        **{f"original_map_{i:02}": p for i, p in enumerate(sorted(original_maps.values()))},
    }
    for name, key in [
        ("km5000noMax.Rdata", "km5000_rdata_sha256"),
        ("RFU.R", "rfu_r_sha256"),
        ("trimerMDSfit_small.Rdata", "trimer_rdata_sha256"),
    ]:
        assert_hash(reference / name, config["reference"][key])
        inputs[name] = reference / name
    out = workspace / "results/rp1_14_repair_qc_v1"
    scientific = identity(
        inputs, (run, vector_hash, vector_audit, conflict_count, stratified_sample)
    )
    fingerprint = start(out, scientific)
    if fingerprint is None:
        validate_completed_rfu_run(out / "parity", expected_provenance=config["reference"])
        return
    frame = pd.read_parquet(prepared / "reused_receptor_assignments.parquet")
    registry = pd.read_csv(prepared / "source_sample_crosswalk.tsv", sep="\t")
    validate_boundary(frame, registry)
    if set(original_maps) != set(registry.source_sample):
        raise ValueError("Source map files do not exactly match the published registry.")
    code = historical.read_text()
    if "cdr3[1:n_map]" not in code or "dd  <- EncodeRepertoire(ff)" not in code:
        raise ValueError("Recovered alignment mechanism not found in historical source.")
    audits, before, repaired_frames = [], [], []
    for row in registry.to_dict("records"):
        saved = pd.read_csv(original_maps[row["source_sample"]], sep="\t")
        selected = frame[frame.sample_id.eq(row["sample_id"])].reset_index(drop=True).copy()
        audit = vector_audit(saved, selected)
        changed = ~saved.cdr3_aa.eq(selected.aminoAcid)
        selected["repair_status"] = np.where(changed, "AA_label_repaired", "AA_label_unchanged")
        audits.append(
            {
                "sample_id": row["sample_id"],
                "donor": row["donor"],
                "visit": row["visit"],
                "compartment": row["compartment"],
                **audit,
                "changed_AA_labels": int(changed.sum()),
                "unchanged_AA_labels": int((~changed).sum()),
                "missing_historical_RFU": int(saved.rfu.isna().sum()),
            }
        )
        before.append(saved)
        repaired_frames.append(selected)
    frame = pd.concat(repaired_frames, ignore_index=True)
    old = pd.concat(before, ignore_index=True)
    conflicts = pd.DataFrame(
        [
            {
                "value": value,
                "before_conflicting_AA_groups": conflict_count(old, "cdr3_aa", source),
                "after_conflicting_AA_groups": conflict_count(frame, "aminoAcid", current),
            }
            for value, source, current in [("RFU", "rfu", "rfu"), ("score", "max_cor", "max_cor")]
        ]
    )
    save(pd.DataFrame(audits), out, "alignment_vector_audit.tsv")
    save(conflicts, out, "repeated_receptor_conflicts.tsv")
    sample, strata, design = stratified_sample(frame)
    save(sample, out, "parity_selected_receptors.tsv")
    save(strata, out, "parity_sampling_strata.tsv")
    json_save(design, out / "parity_design.json")
    assay = out / "parity"
    assay.mkdir(exist_ok=True)
    receptors = pd.DataFrame(
        {
            "input_row_id": [f"parity_{i:04}" for i in range(len(sample))],
            "cell_id": [f"parity_{i:04}" for i in range(len(sample))],
            "chain": "TRB",
            "cdr3aa": sample.aminoAcid,
            "v_call": sample.v_call,
            "source_adapter": "rp1_14_stratified_bounded_QC",
        }
    )
    # Core chunk cache checks input/reference/code/threshold identities on resume.
    result = assign_rfu(
        receptors,
        rfu_dir=reference,
        mode="standard",
        threshold=0.6,
        deduplicate=True,
        chunk_size=400,
        max_workers=1,
        resume=True,
        workdir=assay / "backend",
        wrapper_r_path=REPO / "r/run_rfu_repo.R",
        rscript_bin="/usr/bin/Rscript",
    )
    for name, table in [
        ("receptors", receptors),
        ("unique_sequence_map", result.mapping),
        ("rfu_results_per_sequence", result.per_sequence),
        ("rfu_results_per_row", result.per_row),
    ]:
        table.to_csv(assay / f"{name}.tsv.gz", sep="\t", index=False)
    json_save(
        {**result.provenance, "workflow": "rp1_14_stratified_repair_QC"},
        assay / "run_manifest.json",
    )
    json_save(
        validate_completed_rfu_run(assay, expected_provenance=config["reference"]),
        assay / "validation.json",
    )
    compare = sample.assign(input_row_id=receptors.input_row_id).merge(
        result.per_row[["input_row_id", "rfu_label_nearest", "rfu_score", "rfu_pass_threshold"]],
        on="input_row_id",
        validate="one_to_one",
    )
    compare["same_label"] = compare.rfu.eq(compare.rfu_label_nearest)
    compare["score_error"] = (compare.max_cor - compare.rfu_score).abs()
    compare["same_threshold_status"] = compare.pass_threshold.eq(compare.rfu_pass_threshold)
    save(compare, out, "stratified_reference_parity.tsv")
    summary = []
    for factor in ["all", *STRATA]:
        groups = [("all", compare)] if factor == "all" else compare.groupby(factor, observed=True)
        for name, subset in groups:
            summary.append(
                {
                    "factor": factor,
                    "level": str(name),
                    "n_unique_receptors": len(subset),
                    "label_mismatches": int((~subset.same_label).sum()),
                    "threshold_mismatches": int((~subset.same_threshold_status).sum()),
                    "score_mismatches_at_1e_12": int(subset.score_error.gt(1e-12).sum()),
                    "max_abs_score_difference": float(subset.score_error.max()),
                }
            )
    save(pd.DataFrame(summary), out, "parity_summary.tsv")
    evidence = {
        "source_rows": len(frame),
        "AA_labels_repaired": int(sum(a["changed_AA_labels"] for a in audits)),
        "AA_labels_unchanged": int(sum(a["unchanged_AA_labels"] for a in audits)),
        "RFU_and_score_vectors_identical": True,
        "n_assayed_unique_receptors": len(compare),
        "label_mismatches": int((~compare.same_label).sum()),
        "threshold_mismatches": int((~compare.same_threshold_status).sum()),
        "max_abs_score_difference": float(compare.score_error.max()),
        "score_comparison_tolerance": 1e-12,
        "limitation": "bounded computational parity; does not recover historical invocation/reference execution manifest",
    }
    json_save(evidence, out / "evidence_counts.json")
    save(
        pd.DataFrame(
            [
                {
                    "file": p.name,
                    "sha256": file_sha256(p),
                    "role": "supplementary alignment/parity QC; no biological endpoint rerun",
                }
                for p in sorted(out.glob("*.tsv"))
            ]
        ),
        out,
        "source_table_index.tsv",
    )
    # Include serialized assay outputs in the final manifest, not just the summary.
    outputs = {
        p.name: file_sha256(p) for p in out.iterdir() if p.is_file() and p.name != "completion.json"
    }
    outputs.update(
        {str(p.relative_to(out)): file_sha256(p) for p in assay.iterdir() if p.is_file()}
    )
    status = {"status": "complete", "fingerprint": fingerprint, "outputs": outputs}
    json_save(status, out / "completion.json")
    print(json.dumps(evidence, sort_keys=True), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--reference", type=Path, default=Path("/home/liuyuchen/ext/RFU-official"))
    parser.add_argument(
        "--historical-code",
        type=Path,
        default=Path("/home/liuyuchen/RFU_manuscript_2026/RFU-main/RFU.R"),
    )
    args = parser.parse_args()
    run(args.workspace, args.reference, args.historical_code)


if __name__ == "__main__":
    main()
