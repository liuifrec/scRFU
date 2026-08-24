#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import resource
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from scrfu import __version__, tl
from scrfu.adapters import prepare_receptors
from scrfu.io import file_sha256


def _read_soft_samples(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        text = handle.read()
    rows: list[dict[str, str]] = []
    for block in text.split("^SAMPLE = ")[1:]:
        lines = block.splitlines()
        accession = lines[0].strip()
        values: dict[str, str] = {}
        supplementary = ""
        for line in lines[1:]:
            if line.startswith("!Sample_title = "):
                values["title"] = line.split(" = ", 1)[1]
            elif line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, item = value.split(": ", 1)
                    values[key] = item
            elif "!Sample_supplementary_file_" in line and " = " in line:
                value = line.split(" = ", 1)[1]
                if value.endswith("filtered_contig_annotations.csv.gz"):
                    supplementary = Path(value).name
        if supplementary:
            title = values.get("title", "")
            donor = values.get("subject id") or title.rsplit("_", 1)[0]
            rows.append(
                {
                    "sample_id": accession,
                    "donor_id": donor,
                    "age_group": values.get("age group", ""),
                    "source_assay": values.get("assay", ""),
                    "source_title": title,
                    "filename": supplementary,
                }
            )
    result = pd.DataFrame(rows).sort_values("sample_id", kind="stable").reset_index(drop=True)
    if len(result) != 17 or result[["sample_id", "donor_id", "age_group"]].eq("").any().any():
        raise ValueError("GSE157007 metadata did not resolve to 17 complete VDJ sample mappings.")
    if result["sample_id"].duplicated().any() or result["donor_id"].duplicated().any():
        raise ValueError("Held-out sample and donor mappings must be one-to-one.")
    return result


def _verify_inputs(input_dir: Path, soft: Path, manifest_path: Path) -> tuple[dict[str, Any], str]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected = manifest["hashes_after_download"]
    observed: dict[str, str] = {soft.name: file_sha256(soft)}
    for name in sorted(item for item in expected if item != soft.name):
        observed[name] = file_sha256(input_dir / name)
    mismatches = {
        name: {"expected": expected.get(name), "observed": digest}
        for name, digest in observed.items()
        if expected.get(name) != digest
    }
    if mismatches or set(observed) != set(expected):
        raise ValueError(f"Held-out input hash validation failed: {mismatches}")
    receptor_lines = [
        f"{observed[name]}  {name}\n" for name in sorted(observed) if name != soft.name
    ]
    canonical_hash = hashlib.sha256("".join(receptor_lines).encode()).hexdigest()
    if canonical_hash != manifest["canonical_receptor_manifest_sha256"]:
        raise ValueError("Canonical held-out receptor-manifest hash does not match.")
    return manifest, canonical_hash


def _prepare_receptors(input_dir: Path, samples: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    qc_rows: list[dict[str, Any]] = []
    for sample in samples.itertuples(index=False):
        adapted = prepare_receptors(
            input_dir / sample.filename,
            adapter="cellranger_vdj",
            chain="TRB",
            productive_only=True,
            primary_chain=True,
        )
        frame = adapted.receptors.copy()
        frame["cell_id"] = sample.sample_id + ":" + frame["cell_id"].astype(str)
        frame["input_row_id"] = sample.sample_id + ":" + frame["input_row_id"].astype(str)
        if "clonotype_id" in frame:
            present = frame["clonotype_id"].notna()
            frame.loc[present, "clonotype_id"] = (
                sample.sample_id + ":" + frame.loc[present, "clonotype_id"].astype(str)
            )
        frame["sample_id"] = sample.sample_id
        frame["donor_id"] = sample.donor_id
        frame["age_group"] = sample.age_group
        frames.append(frame)
        qc_rows.append({"sample_id": sample.sample_id, "donor_id": sample.donor_id, **adapted.qc})
    receptors = pd.concat(frames, ignore_index=True)
    if receptors["cell_id"].duplicated().any() or receptors["input_row_id"].duplicated().any():
        raise ValueError("Prepared held-out receptor identifiers are not unique.")
    return receptors, pd.DataFrame(qc_rows)


def _cosine(left: np.ndarray, right: np.ndarray) -> float:
    denominator = np.linalg.norm(left) * np.linalg.norm(right)
    return float(np.dot(left, right) / denominator) if denominator else np.nan


def _sample_structure(
    matrix: pd.DataFrame, samples: pd.DataFrame, representation: str
) -> pd.DataFrame:
    metadata = samples.set_index("sample_id")
    rows: list[dict[str, Any]] = []
    for left_index, left in enumerate(matrix.index):
        for right in matrix.index[left_index + 1 :]:
            rows.append(
                {
                    "representation": representation,
                    "sample_a": left,
                    "sample_b": right,
                    "age_relation": (
                        "within_age_group"
                        if metadata.at[left, "age_group"] == metadata.at[right, "age_group"]
                        else "between_age_group"
                    ),
                    "cosine_similarity": _cosine(
                        matrix.loc[left].to_numpy(dtype=float),
                        matrix.loc[right].to_numpy(dtype=float),
                    ),
                }
            )
    return pd.DataFrame(rows)


def _run_comparators(
    assigned: pd.DataFrame, samples: pd.DataFrame, outdir: Path
) -> tuple[list[dict[str, Any]], list[pd.DataFrame], dict[str, pd.DataFrame]]:
    status: list[dict[str, Any]] = []
    structure: list[pd.DataFrame] = []
    matrices: dict[str, pd.DataFrame] = {}
    methods = ["exact_cdr3", "clonotype", "v_gene", "j_gene", "cdr3_length", "shannon", "simpson"]
    for method in methods:
        try:
            representation = tl.repertoire_representation(
                assigned, sample_key="sample_id", method=method
            )
        except ValueError as exc:
            status.append({"method": method, "status": "skipped", "reason": str(exc)})
            continue
        representation.matrix.to_csv(outdir / f"comparator_{method}.tsv", sep="\t")
        matrices[method] = representation.matrix
        status.append(
            {
                "method": method,
                "status": "completed",
                "samples": len(representation.matrix),
                "features": representation.matrix.shape[1],
            }
        )
        structure.append(_sample_structure(representation.matrix, samples, method))
    unique_cdr3 = assigned["cdr3aa"].nunique()
    status.append(
        {
            "method": "edit_distance",
            "status": "skipped",
            "reason": f"{unique_cdr3} unique CDR3 sequences exceeds frozen quadratic limit 2000",
        }
    )
    return status, structure, matrices


def _deterministic_subsample(frame: pd.DataFrame, fraction: float, seed: int) -> pd.DataFrame:
    selected: list[pd.DataFrame] = []
    for _, sample in frame.groupby("sample_id", observed=True, sort=True):
        count = max(1, int(round(len(sample) * fraction)))
        rank = sample["cell_id"].map(
            lambda value: hashlib.sha256(f"{seed}:{value}".encode()).hexdigest()
        )
        selected.append(sample.assign(_rank=rank).sort_values("_rank", kind="stable").iloc[:count])
    return pd.concat(selected, ignore_index=True).drop(columns="_rank")


def _aligned_cosines(full: pd.DataFrame, perturbed: pd.DataFrame) -> list[float]:
    columns = full.columns.union(perturbed.columns)
    candidate = perturbed.reindex(index=full.index, columns=columns, fill_value=0)
    reference = full.reindex(columns=columns, fill_value=0)
    return [
        _cosine(
            reference.loc[sample].to_numpy(dtype=float), candidate.loc[sample].to_numpy(dtype=float)
        )
        for sample in reference.index
    ]


def _run_stability(
    assigned: pd.DataFrame, full_matrices: dict[str, pd.DataFrame], outdir: Path
) -> None:
    rows: list[dict[str, Any]] = []
    for fraction in (0.5, 0.75):
        for seed in (20260824, 20260825, 20260826):
            subset = _deterministic_subsample(assigned, fraction, seed)
            candidate_matrices: dict[str, pd.DataFrame] = {}
            for policy in ("nearest", "threshold_pass"):
                candidate_matrices[f"rfu_{policy}"] = tl.rfu_pseudobulk(
                    subset,
                    sample_key="sample_id",
                    assignment_policy=policy,
                    weighting="cell",
                    normalize="count",
                ).matrix
            for method in (
                "exact_cdr3",
                "clonotype",
                "v_gene",
                "j_gene",
                "cdr3_length",
                "shannon",
                "simpson",
            ):
                candidate_matrices[method] = tl.repertoire_representation(
                    subset, sample_key="sample_id", method=method
                ).matrix
            for method, full in full_matrices.items():
                for sample, cosine in zip(
                    full.index, _aligned_cosines(full, candidate_matrices[method]), strict=True
                ):
                    rows.append(
                        {
                            "representation": method,
                            "fraction": fraction,
                            "random_seed": seed,
                            "sample_id": sample,
                            "cosine_similarity": cosine,
                        }
                    )
    stability = pd.DataFrame(rows)
    stability.to_csv(outdir / "subsampling_stability.tsv", sep="\t", index=False)
    stability.groupby(["representation", "fraction"], observed=True)["cosine_similarity"].agg(
        ["count", "mean", "median", "std", "min"]
    ).reset_index().to_csv(outdir / "subsampling_stability_summary.tsv", sep="\t", index=False)


def run(args: argparse.Namespace) -> None:
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    manifest, input_hash = _verify_inputs(args.input_dir, args.family_soft, args.manifest)
    samples = _read_soft_samples(args.family_soft)
    receptors, sample_qc = _prepare_receptors(args.input_dir, samples)
    reference = tl.validate_frozen_reference(json.loads(args.frozen_reference.read_text()))

    started_cpu = time.process_time()
    started_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.perf_counter()
    wrapper = Path(__file__).resolve().parents[1] / "r" / "run_rfu_repo.R"
    result = tl.call_rfu_table(
        receptors,
        rfu_dir=args.rfu_dir,
        mode="standard",
        threshold=0.6,
        deduplicate=True,
        chunk_size=5000,
        max_workers=2,
        executor="process",
        resume=True,
        force_recompute=False,
        workdir=outdir / "rfu_cache",
        wrapper_r_path=wrapper,
    )
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    timing = {
        "wall_time_seconds": time.perf_counter() - started,
        "python_cpu_time_seconds": time.process_time() - started_cpu,
        "child_cpu_time_seconds": (children.ru_utime + children.ru_stime)
        - (started_children.ru_utime + started_children.ru_stime),
        "python_peak_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "backend_peak_rss_kb": result.provenance.get("process_max_rss_kb"),
    }
    assigned = result.per_row
    assigned.to_csv(outdir / "assigned_rows.tsv.gz", sep="\t", index=False)
    samples.to_csv(outdir / "sample_metadata.tsv", sep="\t", index=False)
    sample_qc.to_csv(outdir / "sample_receptor_qc.tsv", sep="\t", index=False)

    structure_tables: list[pd.DataFrame] = []
    transfer_summary: list[dict[str, Any]] = []
    full_matrices: dict[str, pd.DataFrame] = {}
    for policy in ("nearest", "threshold_pass"):
        transfer = tl.transfer_cohort(
            assigned,
            reference,
            cohort_label="GSE157007_heldout_aging_frailty",
            sample_key="sample_id",
            assignment_policy=policy,
            weighting="cell",
            groupby="age_group",
            observed_reference_id=reference.immutable_reference_id,
        )
        transfer.coverage.to_csv(outdir / f"reference_coverage_{policy}.tsv", sep="\t", index=False)
        transfer.score_distribution.to_csv(
            outdir / f"score_distribution_{policy}.tsv", sep="\t", index=False
        )
        transfer.rfu_summary.to_csv(outdir / f"rfu_summary_{policy}.tsv", sep="\t", index=False)
        transfer.sample_matrix.to_csv(outdir / f"rfu_pseudobulk_{policy}.tsv", sep="\t")
        full_matrices[f"rfu_{policy}"] = transfer.sample_matrix
        structure_tables.append(_sample_structure(transfer.sample_matrix, samples, f"rfu_{policy}"))
        transfer_summary.append(
            {
                "assignment_policy": policy,
                "receptor_rows": len(assigned),
                "unique_cdr3": assigned["cdr3aa"].nunique(),
                "samples": len(transfer.sample_matrix),
                "rfus": transfer.sample_matrix.shape[1],
                "threshold_pass_fraction": float(assigned["pass_thr"].astype("boolean").mean()),
            }
        )

    comparator_status, comparator_structure, comparator_matrices = _run_comparators(
        assigned, samples, outdir
    )
    full_matrices.update(comparator_matrices)
    _run_stability(assigned, full_matrices, outdir)
    structure = pd.concat([*structure_tables, *comparator_structure], ignore_index=True)
    structure.to_csv(outdir / "sample_structure.tsv", sep="\t", index=False)
    structure.groupby(["representation", "age_relation"], observed=True)["cosine_similarity"].agg(
        ["count", "mean", "median", "std"]
    ).reset_index().to_csv(outdir / "sample_structure_summary.tsv", sep="\t", index=False)
    pd.DataFrame(comparator_status).to_csv(outdir / "comparator_status.tsv", sep="\t", index=False)
    pd.DataFrame(transfer_summary).to_csv(outdir / "transfer_summary.tsv", sep="\t", index=False)

    run_manifest = {
        "dataset_label": manifest["dataset_label"],
        "input_manifest_sha256": input_hash,
        "frozen_reference_id": reference.immutable_reference_id,
        "software_version": __version__,
        "parameters": {
            "chain": "TRB",
            "productive_only": True,
            "primary_chain": True,
            "threshold": 0.6,
            "chunk_size": 5000,
            "workers": 2,
            "assignment_policies": ["nearest", "threshold_pass"],
        },
        "source_rows": int(sample_qc["source_row_count"].sum()),
        "receptor_rows": len(receptors),
        "unique_cdr3": receptors["cdr3aa"].nunique(),
        "sample_count": len(samples),
        "donor_count": samples["donor_id"].nunique(),
        "age_group_count": samples["age_group"].nunique(),
        "rfu_provenance": result.provenance,
        "timing": timing,
        "scope": "Prespecified held-out technical transfer; age-group summaries are descriptive.",
    }
    (outdir / "run_manifest.json").write_text(
        json.dumps(run_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(run_manifest, indent=2, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run preregistered GSE157007 held-out transfer.")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--family-soft", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--frozen-reference", type=Path, required=True)
    parser.add_argument("--rfu-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser


if __name__ == "__main__":
    run(build_parser().parse_args())
