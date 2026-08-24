#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import resource
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from scrfu import __version__, tl
from scrfu.adapters import adapt_wells_tcr_ir
from scrfu.io import file_sha256, read_h5ad_obs, write_receptor_cache
from scrfu.wells import read_wells_receptors_h5ad, source_fingerprint

METADATA_COLUMNS = [
    "donor_id",
    "donor_age",
    "sex",
    "cmv",
    "disease",
    "tissue",
    "cell_type",
    "cell_state",
    "library_id",
]


def _json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def _selected_hash(values: pd.Series | pd.Index) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(str(value).encode())
        digest.update(b"\n")
    return digest.hexdigest()


def _directory_size(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def _proportional_stratified_sample(
    obs: pd.DataFrame, *, n: int, strata: list[str], seed: int
) -> pd.DataFrame:
    grouped = obs.groupby(strata, observed=True, dropna=False, sort=True)
    sizes = grouped.size().rename("size").reset_index()
    exact = sizes["size"] * n / len(obs)
    sizes["quota"] = np.floor(exact).astype(int)
    sizes["remainder"] = exact - sizes["quota"]
    remaining = n - int(sizes["quota"].sum())
    sizes = sizes.sort_values(["remainder", *strata], ascending=[False, *([True] * len(strata))])
    sizes.loc[sizes.index[:remaining], "quota"] += 1
    quotas = {
        tuple(row[column] for column in strata): int(row["quota"]) for _, row in sizes.iterrows()
    }
    selected: list[pd.DataFrame] = []
    for key, frame in grouped:
        key_tuple = key if isinstance(key, tuple) else (key,)
        quota = quotas[key_tuple]
        if quota:
            ranked = frame.assign(
                _rank=frame["cell_id"].map(
                    lambda cell: hashlib.sha256(f"{seed}:{cell}".encode()).hexdigest()
                )
            ).sort_values("_rank", kind="stable")
            selected.append(ranked.iloc[:quota].drop(columns="_rank"))
    result = pd.concat(selected, ignore_index=False).sort_index(kind="stable")
    if len(result) != n:
        raise RuntimeError(f"Stratified sampler selected {len(result)} rows; expected {n}.")
    return result


def _sampling_qc(frame: pd.DataFrame, *, label: str, seed: int, source_rows: int) -> dict[str, Any]:
    return {
        "dataset_label": label,
        "strategy": "deterministic_random_without_replacement",
        "random_seed": seed,
        "selected_cell_count": len(frame),
        "selected_cell_sha256": _selected_hash(frame["cell_id"]),
        "source_cell_count": source_rows,
        "donor_count": int(frame["donor_id"].nunique()),
        "cell_type_count": int(frame["cell_type"].nunique()),
        "tissue_count": int(frame["tissue"].nunique()),
        "library_count": int(frame["library_id"].nunique()),
    }


def _stratum_total_variation(
    source: pd.DataFrame, sample: pd.DataFrame, strata: list[str]
) -> float:
    expected = source.groupby(strata, observed=True, dropna=False).size() / len(source)
    observed = sample.groupby(strata, observed=True, dropna=False).size() / len(sample)
    index = expected.index.union(observed.index)
    return float(
        0.5
        * np.abs(
            expected.reindex(index, fill_value=0) - observed.reindex(index, fill_value=0)
        ).sum()
    )


def _measure_call(function: Any, *args: Any, **kwargs: Any) -> tuple[Any, dict[str, float]]:
    before_cpu = time.process_time()
    before_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.perf_counter()
    result = function(*args, **kwargs)
    after_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return result, {
        "wall_time_seconds": time.perf_counter() - started,
        "python_cpu_time_seconds": time.process_time() - before_cpu,
        "child_cpu_time_seconds": (after_children.ru_utime + after_children.ru_stime)
        - (before_children.ru_utime + before_children.ru_stime),
        "peak_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }


def _assignment_hash(frame: pd.DataFrame) -> str:
    columns = ["input_row_id", "rfu_id", "rfu_label", "rfu_score", "pass_thr"]
    normalized = frame.loc[:, columns].sort_values("input_row_id", kind="stable")
    return hashlib.sha256(normalized.to_csv(index=False).encode()).hexdigest()


def _assert_equal(reference: pd.DataFrame, candidate: pd.DataFrame) -> None:
    columns = [
        "rfu_id",
        "rfu_label",
        "rfu_score",
        "pass_thr",
        "eligibility_status",
        "reference_coverage_status",
    ]
    left = reference.set_index("input_row_id").sort_index()
    right = candidate.set_index("input_row_id").reindex(left.index)
    if right.index.has_duplicates or right[columns].isna().all(axis=1).any():
        raise AssertionError("Candidate RFU output has missing or duplicate input rows.")
    for column in (
        "rfu_id",
        "rfu_label",
        "pass_thr",
        "eligibility_status",
        "reference_coverage_status",
    ):
        left_values = left[column].astype("string").fillna("<NA>")
        right_values = right[column].astype("string").fillna("<NA>")
        if not left_values.equals(right_values):
            raise AssertionError(f"RFU invariance failed for {column}.")
    if not np.allclose(
        pd.to_numeric(left["rfu_score"], errors="coerce"),
        pd.to_numeric(right["rfu_score"], errors="coerce"),
        rtol=0,
        atol=1e-12,
        equal_nan=True,
    ):
        raise AssertionError("RFU invariance failed for rfu_score.")


def _run_rfu(
    receptors: pd.DataFrame,
    *,
    rfu_dir: Path,
    wrapper: Path,
    workdir: Path,
    chunk_size: int | None,
    workers: int,
    resume: bool,
    force: bool,
) -> tuple[Any, dict[str, float]]:
    return _measure_call(
        tl.call_rfu_table,
        receptors,
        rfu_dir=rfu_dir,
        mode="standard",
        threshold=0.6,
        deduplicate=True,
        chunk_size=chunk_size,
        max_workers=workers,
        executor="process",
        resume=resume,
        force_recompute=force,
        workdir=workdir,
        wrapper_r_path=wrapper,
    )


def _benchmark_row(
    *,
    label: str,
    phase: str,
    result: Any,
    timing: dict[str, float],
    receptors: pd.DataFrame,
    workdir: Path,
    chunk_size: int | None,
    workers: int,
    input_hash: str,
    rfu_reference: str,
) -> dict[str, Any]:
    per_row = result.per_row
    eligible = per_row["eligibility_status"].eq("eligible")
    provenance = result.provenance
    return {
        "dataset_label": label,
        "phase": phase,
        "input_hash": input_hash,
        "rfu_reference_identifier": rfu_reference,
        "software_version": __version__,
        "parameters": json.dumps(
            {"threshold": 0.6, "chunk_size": chunk_size, "workers": workers, "executor": "process"},
            sort_keys=True,
        ),
        "random_seed": 20260824,
        "source_cells": int(label.removeprefix("wells_")),
        "receptor_rows": len(receptors),
        "eligible_rows": int(eligible.sum()),
        "unique_cdr3_queries": int(provenance["unique_query_count"]),
        "deduplication_ratio": provenance["deduplication_ratio"],
        "rfus_observed": int(per_row.loc[eligible, "rfu_label"].nunique()),
        "threshold_pass_rate": float(per_row.loc[eligible, "pass_thr"].astype("boolean").mean()),
        "chunk_count": int(provenance["chunk_count"]),
        "worker_count": workers,
        "executor": "process",
        **timing,
        "backend_peak_rss_kb": provenance.get("process_max_rss_kb"),
        "reused_chunk_count": int(provenance["reused_chunk_count"]),
        "recomputed_chunk_count": int(provenance["recomputed_chunk_count"]),
        "output_size_bytes": _directory_size(workdir),
        "assignment_sha256": _assignment_hash(per_row),
    }


def _representation(frame: pd.DataFrame, policy: str) -> pd.DataFrame:
    return tl.rfu_pseudobulk(
        frame,
        sample_key="library_id",
        assignment_policy=policy,
        weighting="cell",
        normalize="proportion",
    ).matrix


def _stability_rows(
    reference: pd.DataFrame,
    perturbed: pd.DataFrame,
    *,
    perturbation: str,
    fraction: float,
    seed: int,
    policy: str,
    coverage_change: float,
) -> pd.DataFrame:
    metrics = tl.benchmark_representation_stability(reference, perturbed, top_k=10).metrics
    metrics["perturbation"] = perturbation
    metrics["fraction"] = fraction
    metrics["random_seed"] = seed
    metrics["assignment_policy"] = policy
    metrics["reference_coverage_change"] = coverage_change
    return metrics


def _phenotype_coupling_stability(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for policy in ("nearest", "threshold_pass"):
        reference = tl.rfu_phenotype_coupling(
            frame,
            phenotype_key="cell_type",
            assignment_policy=policy,
            weighting="cell",
        )
        reference_vector = reference.set_index(["rfu_label", "cell_type"])[
            "phenotype_specific_proportion"
        ]
        reference_dominant = reference.drop_duplicates("rfu_label").set_index("rfu_label")[
            "dominant_phenotype"
        ]
        for unit in ("cell", "sequence"):
            for fraction in (0.25, 0.5, 0.75, 1.0):
                seeds = (20260824,) if fraction == 1.0 else (17, 29, 43)
                for seed in seeds:
                    subset = tl.deterministic_subsample(
                        frame,
                        unit=unit,
                        fraction=fraction,
                        random_state=seed,
                        cell_col="cell_id",
                        sequence_col="cdr3aa",
                    )
                    perturbed = tl.rfu_phenotype_coupling(
                        subset,
                        phenotype_key="cell_type",
                        assignment_policy=policy,
                        weighting="cell",
                    )
                    perturbed_vector = perturbed.set_index(["rfu_label", "cell_type"])[
                        "phenotype_specific_proportion"
                    ]
                    index = reference_vector.index.union(perturbed_vector.index)
                    left = reference_vector.reindex(index, fill_value=0.0)
                    right = perturbed_vector.reindex(index, fill_value=0.0)
                    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
                    perturbed_dominant = perturbed.drop_duplicates("rfu_label").set_index(
                        "rfu_label"
                    )["dominant_phenotype"]
                    shared = reference_dominant.index.intersection(perturbed_dominant.index)
                    rows.append(
                        {
                            "assignment_policy": policy,
                            "perturbation": f"{unit}_subsampling",
                            "fraction": fraction,
                            "random_seed": seed,
                            "spearman": float(left.corr(right, method="spearman")),
                            "cosine": float(np.dot(left, right) / denominator)
                            if denominator
                            else np.nan,
                            "mean_absolute_error": float(np.abs(left - right).mean()),
                            "shared_rfu_count": len(shared),
                            "dominant_phenotype_agreement": float(
                                reference_dominant.loc[shared]
                                .eq(perturbed_dominant.loc[shared])
                                .mean()
                            )
                            if len(shared)
                            else np.nan,
                        }
                    )
    return pd.DataFrame(rows)


def _run_robustness(
    frame: pd.DataFrame,
    outdir: Path,
    *,
    input_hash: str,
    rfu_reference: str,
) -> None:
    rows: list[pd.DataFrame] = []
    for policy in ("nearest", "threshold_pass"):
        reference = _representation(frame, policy)
        full_coverage = float(frame["pass_thr"].astype("boolean").mean())
        for unit in ("cell", "sequence"):
            for fraction in (0.25, 0.5, 0.75, 1.0):
                seeds = (20260824,) if fraction == 1.0 else (17, 29, 43)
                for seed in seeds:
                    subset = tl.deterministic_subsample(
                        frame,
                        unit=unit,
                        fraction=fraction,
                        random_state=seed,
                        cell_col="cell_id",
                        sequence_col="cdr3aa",
                    )
                    rows.append(
                        _stability_rows(
                            reference,
                            _representation(subset, policy),
                            perturbation=f"{unit}_subsampling",
                            fraction=fraction,
                            seed=seed,
                            policy=policy,
                            coverage_change=float(subset["pass_thr"].astype("boolean").mean())
                            - full_coverage,
                        )
                    )
        counts = tl.rfu_pseudobulk(
            frame,
            sample_key="library_id",
            assignment_policy=policy,
            weighting="cell",
            normalize="count",
        ).matrix
        proportions = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
        for fraction in (0.25, 0.5, 0.75, 1.0):
            seeds = (20260824,) if fraction == 1.0 else (17, 29, 43)
            depths = np.ceil(counts.sum(axis=1) * fraction).astype(int)
            for seed in seeds:
                resampled = tl.multinomial_abundance_resample(
                    counts, depth=depths, random_state=seed
                )
                perturbed = resampled.div(resampled.sum(axis=1).replace(0, np.nan), axis=0).fillna(
                    0
                )
                rows.append(
                    _stability_rows(
                        proportions,
                        perturbed,
                        perturbation="multinomial_abundance_resampling",
                        fraction=fraction,
                        seed=seed,
                        policy=policy,
                        coverage_change=0.0,
                    )
                )
    robustness = pd.concat(rows, ignore_index=True)
    provenance = {
        "dataset_label": "wells_25000",
        "input_hash": input_hash,
        "rfu_reference_identifier": rfu_reference,
        "software_version": __version__,
        "parameters": json.dumps(
            {
                "fractions": [0.25, 0.5, 0.75, 1.0],
                "seeds": [17, 29, 43],
                "top_k": 10,
                "multinomial_abundance_is_physical_read_downsampling": False,
            },
            sort_keys=True,
        ),
    }
    for column, value in provenance.items():
        robustness[column] = value
    robustness.to_csv(outdir / "figure1_robustness_source.tsv", sep="\t", index=False)
    coupling = _phenotype_coupling_stability(frame)
    for column, value in provenance.items():
        coupling[column] = value
    coupling.to_csv(outdir / "phenotype_coupling_stability.tsv", sep="\t", index=False)
    threshold = tl.threshold_sensitivity(frame, (0.5, 0.55, 0.6, 0.65, 0.7), groupby="library_id")
    for column, value in provenance.items():
        threshold[column] = value
    threshold.to_csv(outdir / "threshold_sensitivity.tsv", sep="\t", index=False)


def _run_single_cell_analysis(frame: pd.DataFrame, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    policy_summary: list[dict[str, Any]] = []
    for policy in ("nearest", "threshold_pass"):
        policy_dir = outdir / policy
        policy_dir.mkdir(parents=True, exist_ok=True)
        metrics = tl.rfu_metrics(
            frame,
            groupby="cell_type",
            weighting="cell",
            assignment_policy=policy,
            donor_col="donor_id",
            sample_col="library_id",
        )
        metrics.to_csv(policy_dir / "rfu_metrics_by_cell_type.tsv", sep="\t", index=False)
        repertoire = tl.repertoire_metrics(frame, groupby="library_id", weighting="cell")
        repertoire.to_csv(
            policy_dir / "conventional_repertoire_by_library.tsv", sep="\t", index=False
        )
        coupling = tl.rfu_phenotype_coupling(
            frame,
            phenotype_key="cell_type",
            sample_key="library_id",
            assignment_policy=policy,
            weighting="cell",
        )
        coupling.to_csv(policy_dir / "rfu_phenotype_coupling.tsv", sep="\t", index=False)
        for sample_key in ("library_id", "donor_id"):
            pseudobulk = tl.rfu_pseudobulk(
                frame,
                sample_key=sample_key,
                assignment_policy=policy,
                weighting="cell",
                normalize="count",
            )
            pseudobulk.matrix.to_csv(policy_dir / f"rfu_pseudobulk_by_{sample_key}.tsv", sep="\t")
            for metric in ("jaccard", "cosine"):
                tl.rfu_overlap(pseudobulk, metric=metric).matrix.to_csv(
                    policy_dir / f"rfu_{metric}_by_{sample_key}.tsv", sep="\t"
                )
        policy_summary.append(
            {
                "assignment_policy": policy,
                "assigned_cells": int(
                    frame["rfu_label"].notna().sum()
                    if policy == "nearest"
                    else frame["pass_thr"].astype("boolean").fillna(False).sum()
                ),
                "observed_rfus": int(
                    frame.loc[
                        frame["rfu_label"].notna()
                        if policy == "nearest"
                        else frame["pass_thr"].astype("boolean").fillna(False),
                        "rfu_label",
                    ].nunique()
                ),
                "phenotype_coupling_rows": len(coupling),
                "cell_types": int(frame["cell_type"].nunique()),
                "donors": int(frame["donor_id"].nunique()),
                "libraries": int(frame["library_id"].nunique()),
            }
        )
    pd.DataFrame(policy_summary).to_csv(
        outdir / "assignment_policy_summary.tsv", sep="\t", index=False
    )
    comparator_status: list[dict[str, Any]] = []
    requirements = {
        "exact_cdr3": "cdr3aa",
        "clonotype": "clonotype_id",
        "v_gene": "v_call",
        "j_gene": "j_call",
        "cdr3_length": "cdr3aa",
        "shannon": "cdr3aa",
        "simpson": "cdr3aa",
    }
    for method, required_column in requirements.items():
        if required_column not in frame:
            comparator_status.append(
                {
                    "method": method,
                    "status": "skipped",
                    "reason": f"input lacks {required_column}",
                }
            )
            continue
        representation = tl.repertoire_representation(frame, sample_key="library_id", method=method)
        representation.matrix.to_csv(outdir / f"comparator_{method}.tsv", sep="\t")
        comparator_status.append(
            {
                "method": method,
                "status": "completed",
                "reason": "",
                "sample_count": representation.matrix.shape[0],
                "feature_count": representation.matrix.shape[1],
            }
        )
    unique_sequences = int(frame["cdr3aa"].nunique())
    if unique_sequences <= 2000:
        representation = tl.repertoire_representation(
            frame, sample_key="library_id", method="edit_distance"
        )
        representation.matrix.to_csv(outdir / "comparator_edit_distance.tsv", sep="\t")
        comparator_status.append(
            {
                "method": "edit_distance",
                "status": "completed",
                "reason": "",
                "sample_count": representation.matrix.shape[0],
                "feature_count": representation.matrix.shape[1],
            }
        )
    else:
        comparator_status.append(
            {
                "method": "edit_distance",
                "status": "skipped",
                "reason": (
                    f"{unique_sequences} unique sequences exceed the bounded quadratic "
                    "comparison limit of 2000"
                ),
            }
        )
    pd.DataFrame(comparator_status).to_csv(outdir / "comparator_status.tsv", sep="\t", index=False)


def run(args: argparse.Namespace) -> None:
    source = args.input.expanduser().resolve()
    rfu_dir = args.rfu_dir.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    wrapper = Path(__file__).resolve().parents[1] / "r" / "run_rfu_repo.R"

    all_obs, metadata_timing = _measure_call(read_h5ad_obs, source, columns=METADATA_COLUMNS)
    all_obs = all_obs.rename_axis("cell_id").reset_index()
    source_state = source_fingerprint(source)
    rfu_reference = file_sha256(rfu_dir / "km5000noMax.Rdata")
    strata = ["donor_id", "tissue", "cell_type"]
    sizes = sorted(set(args.sizes))
    sampling_rows: list[dict[str, Any]] = []
    samples: dict[int, pd.DataFrame] = {}
    for size in sizes:
        sample = tl.deterministic_subsample(
            all_obs,
            unit="cell",
            n=size,
            random_state=args.seed,
            cell_col="cell_id",
            preserve_order=True,
        )
        samples[size] = sample
        qc = _sampling_qc(sample, label=f"wells_{size}", seed=args.seed, source_rows=len(all_obs))
        qc["stratum_total_variation_from_source"] = _stratum_total_variation(
            all_obs, sample, strata
        )
        sampling_rows.append(qc)
    largest = sizes[-1]
    stratified = _proportional_stratified_sample(
        all_obs,
        n=largest,
        strata=strata,
        seed=args.seed,
    )
    stratified_qc = _sampling_qc(
        stratified,
        label=f"wells_{largest}_stratified",
        seed=args.seed,
        source_rows=len(all_obs),
    )
    stratified_qc["strategy"] = "proportional_donor_tissue_cell_type_stratified_hash_rank"
    stratified_qc["stratum_total_variation_from_source"] = _stratum_total_variation(
        all_obs, stratified, strata
    )
    sampling_rows.append(stratified_qc)

    for row in sampling_rows:
        row.update(
            {
                "input_hash": source_state["fingerprint"],
                "rfu_reference_identifier": rfu_reference,
                "software_version": __version__,
                "parameters": json.dumps(
                    {"strata": strata, "preserve_source_order": True}, sort_keys=True
                ),
            }
        )

    benchmark_rows: list[dict[str, Any]] = []
    reference_results: dict[int, Any] = {}
    merged_largest: pd.DataFrame | None = None
    for size in sizes:
        label = f"wells_{size}"
        subset_dir = outdir / "subsets" / label
        selected = samples[size]
        receptor_data, extraction_timing = _measure_call(
            read_wells_receptors_h5ad,
            source,
            obs_columns=METADATA_COLUMNS,
            selected_obs_names=selected["cell_id"].tolist(),
        )
        adapted, adapter_timing = _measure_call(
            adapt_wells_tcr_ir, receptor_data, chain="TRB", primary_chain=True
        )
        cache_dir = subset_dir / "receptor_cache"
        cache_manifest, cache_timing = _measure_call(
            write_receptor_cache,
            cache_dir,
            adapted.receptors,
            receptor_data.obs,
            source_adapter=adapted.adapter_name,
            source_adapter_version=adapted.adapter_version,
            source_format="wells_h5ad",
            source_path=source,
            adapter_qc=adapted.qc,
            adapter_configuration={"chain": "TRB", "primary_chain": True},
            selected_metadata_columns=METADATA_COLUMNS,
            source_atlas_dimensions=receptor_data.atlas_shape,
            force=True,
        )
        input_hash = cache_manifest["files"]["receptors"]["sha256"]
        sample_qc = next(row for row in sampling_rows if row["dataset_label"] == label)
        sample_qc.update(
            {
                "source_dimensions": list(receptor_data.atlas_shape),
                "receptor_bearing_cell_count": int(adapted.receptors["cell_id"].nunique()),
                "productive_trb_count": int(adapted.receptors["productive"].fillna(False).sum()),
                "receptor_cache_size_bytes": _directory_size(cache_dir),
                "targeted_extraction_timing": extraction_timing,
                "adapter_timing": adapter_timing,
                "cache_preparation_timing": cache_timing,
            }
        )
        _json(subset_dir / "sampling_manifest.json", sample_qc)

        baseline_dir = subset_dir / "rfu_unthreaded_unchunked"
        baseline, timing = _run_rfu(
            adapted.receptors,
            rfu_dir=rfu_dir,
            wrapper=wrapper,
            workdir=baseline_dir,
            chunk_size=None,
            workers=1,
            resume=False,
            force=True,
        )
        reference_results[size] = baseline
        baseline.per_row.to_csv(subset_dir / "rfu_results_per_row.tsv.gz", sep="\t", index=False)
        benchmark_rows.append(
            _benchmark_row(
                label=label,
                phase="fresh",
                result=baseline,
                timing=timing,
                receptors=adapted.receptors,
                workdir=baseline_dir,
                chunk_size=None,
                workers=1,
                input_hash=input_hash,
                rfu_reference=rfu_reference,
            )
        )
        unique_queries = int(baseline.provenance["unique_query_count"])
        for chunk_size in args.chunk_sizes:
            if unique_queries < chunk_size:
                continue
            for workers in (1, 2):
                workdir = subset_dir / f"rfu_chunk_{chunk_size}_workers_{workers}"
                fresh, fresh_timing = _run_rfu(
                    adapted.receptors,
                    rfu_dir=rfu_dir,
                    wrapper=wrapper,
                    workdir=workdir,
                    chunk_size=chunk_size,
                    workers=workers,
                    resume=False,
                    force=True,
                )
                _assert_equal(baseline.per_row, fresh.per_row)
                benchmark_rows.append(
                    _benchmark_row(
                        label=label,
                        phase="fresh",
                        result=fresh,
                        timing=fresh_timing,
                        receptors=adapted.receptors,
                        workdir=workdir,
                        chunk_size=chunk_size,
                        workers=workers,
                        input_hash=input_hash,
                        rfu_reference=rfu_reference,
                    )
                )
                resumed, resumed_timing = _run_rfu(
                    adapted.receptors,
                    rfu_dir=rfu_dir,
                    wrapper=wrapper,
                    workdir=workdir,
                    chunk_size=chunk_size,
                    workers=workers,
                    resume=True,
                    force=False,
                )
                _assert_equal(baseline.per_row, resumed.per_row)
                benchmark_rows.append(
                    _benchmark_row(
                        label=label,
                        phase="resumed",
                        result=resumed,
                        timing=resumed_timing,
                        receptors=adapted.receptors,
                        workdir=workdir,
                        chunk_size=chunk_size,
                        workers=workers,
                        input_hash=input_hash,
                        rfu_reference=rfu_reference,
                    )
                )
        if size == largest:
            shuffled = tl.shuffle_input_order(
                adapted.receptors, random_state=args.seed, preserve_index=False
            )
            shuffled_result, shuffled_timing = _run_rfu(
                shuffled,
                rfu_dir=rfu_dir,
                wrapper=wrapper,
                workdir=subset_dir / "rfu_shuffled",
                chunk_size=1000,
                workers=2,
                resume=False,
                force=True,
            )
            _assert_equal(baseline.per_row, shuffled_result.per_row)
            benchmark_rows.append(
                _benchmark_row(
                    label=label,
                    phase="shuffled_order",
                    result=shuffled_result,
                    timing=shuffled_timing,
                    receptors=adapted.receptors,
                    workdir=subset_dir / "rfu_shuffled",
                    chunk_size=1000,
                    workers=2,
                    input_hash=input_hash,
                    rfu_reference=rfu_reference,
                )
            )
            metadata = receptor_data.obs.rename_axis("cell_id").reset_index()
            merged_largest = baseline.per_row.merge(
                metadata, on="cell_id", how="left", validate="many_to_one", sort=False
            )

    pd.DataFrame(sampling_rows).to_csv(outdir / "wells_sampling_summary.tsv", sep="\t", index=False)
    benchmark = pd.DataFrame(benchmark_rows)
    benchmark.to_csv(outdir / "figure1_scaling_source.tsv", sep="\t", index=False)
    if merged_largest is None:
        raise RuntimeError("Largest Wells result was not assembled.")
    _run_robustness(
        merged_largest,
        outdir,
        input_hash=input_hash,
        rfu_reference=rfu_reference,
    )
    _run_single_cell_analysis(merged_largest, outdir / "single_cell_analysis")
    _json(
        outdir / "run_manifest.json",
        {
            "dataset_label": "Wells bounded deterministic validation",
            "source": source_state,
            "rfu_reference_identifier": rfu_reference,
            "rfu_artifacts": {
                name: {
                    "sha256": file_sha256(rfu_dir / name),
                    "size_bytes": (rfu_dir / name).stat().st_size,
                }
                for name in ("RFU.R", "trimerMDSfit_small.Rdata", "km5000noMax.Rdata")
            },
            "software_version": __version__,
            "parameters": {
                "sizes": sizes,
                "chunk_sizes": args.chunk_sizes,
                "random_seed": args.seed,
                "threshold": 0.6,
                "metadata_columns": METADATA_COLUMNS,
            },
            "metadata_read_timing": metadata_timing,
            "source_tables": [
                "wells_sampling_summary.tsv",
                "figure1_scaling_source.tsv",
                "figure1_robustness_source.tsv",
                "threshold_sensitivity.tsv",
                "single_cell_analysis/assignment_policy_summary.tsv",
            ],
            "scientific_scope": (
                "Bounded technical and descriptive single-cell validation; deterministic random "
                "subsets are not asserted to be biologically representative."
            ),
        },
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run bounded Month 2 validation on a Wells H5AD.")
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--rfu-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--sizes", type=int, nargs="+", default=[1000, 10000, 25000])
    parser.add_argument("--chunk-sizes", type=int, nargs="+", default=[1000, 5000, 20000])
    parser.add_argument("--seed", type=int, default=20260824)
    return parser


if __name__ == "__main__":
    run(build_parser().parse_args())
