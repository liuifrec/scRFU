from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pandas as pd

from .benchmark import (
    StabilityBenchmarkResult,
    benchmark_representation_stability,
    deterministic_subsample,
    donor_leave_one_out,
    multinomial_abundance_resample,
    shuffle_input_order,
    threshold_sensitivity,
)
from .comparators import (
    ComparatorRepresentation,
    list_comparators,
    register_comparator,
    repertoire_representation,
)
from .completed_run import validate_completed_rfu_run
from .diagnostics import reference_coverage
from .downstream import (
    RFUOverlapResult,
    RFUPseudobulkResult,
    rfu_overlap,
    rfu_phenotype_coupling,
    rfu_pseudobulk,
)
from .longitudinal import (
    LongitudinalCompartmentResult,
    LongitudinalDesign,
    LongitudinalDynamicsResult,
    LongitudinalResamplingResult,
    RFULongitudinalResult,
    bootstrap_longitudinal_statistic,
    donor_retrieval,
    longitudinal_compartment_comparison,
    longitudinal_similarity,
    permute_longitudinal_labels,
    rfu_longitudinal_dynamics,
    rfu_longitudinal_matrix,
    summarize_longitudinal_similarity,
    validate_longitudinal_design,
)
from .repertoire import repertoire_metrics
from .scverse import (
    ScverseRFUResult,
    SummaryPolicy,
    assign_rfu_native,
    concat_scrfu,
    validate_scrfu_schema,
)
from .sequence import rfu_sequence_matrix
from .summary import aggregate_rfu, rfu_metrics, rfu_summary
from .tl_rfu_repo import RFUTableResult, call_rfu_repo, call_rfu_table
from .transfer import (
    CohortHarmonizationResult,
    FrozenRFUReference,
    HeldOutValidationManifest,
    TransferCohortResult,
    create_heldout_validation_manifest,
    harmonize_cohort_metadata,
    transfer_cohort,
    validate_frozen_reference,
)
from .validation import validate_airr
from .vdjdb import (
    AntigenContextResult,
    AntigenPermutationResult,
    VDJdbEvidenceSummary,
    VDJdbReference,
    annotate_vdjdb,
    compare_antigen_groupings,
    global_antigen_coherence,
    load_vdjdb_reference,
    normalize_vdjdb_cdr3,
    normalize_vdjdb_v_gene,
    rfu_antigen_abundance,
    rfu_antigen_coherence,
    rfu_antigen_permutation_test,
    summarize_antigen_context,
    summarize_vdjdb_evidence,
    validate_vdjdb_reference,
)

PathLike = str | Path


def assign_rfu(
    data: Any,
    *,
    rfu_dir: PathLike | None = None,
    mode: str = "standard",
    threshold: float = 0.6,
    deduplicate: bool = True,
    chunk_size: int | None = None,
    max_workers: int = 1,
    executor: str = "process",
    resume: bool = True,
    force_recompute: bool = False,
    wrapper_r_path: PathLike = "r/run_rfu_repo.R",
    rscript_bin: str = "Rscript",
    extra_r_args: Sequence[str] | None = None,
    workdir: PathLike | None = None,
    inplace: bool = True,
    key_added: str = "scrfu",
    airr_mod: str = "airr",
    airr_key: str = "airr",
    chain_idx_key: str = "chain_indices",
    cell_summary: bool = False,
    summary_policy: SummaryPolicy = "ambiguity_aware",
    summary_key_added: str | None = None,
) -> RFUTableResult | ScverseRFUResult | None:
    """Assign canonical TRB RFUs from a table, AnnData, or MuData object.

    DataFrame inputs delegate directly to :func:`call_rfu_table`. Current
    Scirpy AIRR objects are annotated chain-by-chain in ``obsm[key_added]``;
    MuData results are stored in the selected AIRR modality. Native execution
    never reads ``X`` and preserves every source receptor chain.

    Parameters
    ----------
    data
        Canonical receptor DataFrame, AnnData with ``obsm[airr_key]``, or
        MuData containing ``mod[airr_mod].obsm[airr_key]``.
    inplace
        For AnnData/MuData, populate the native chain/provenance namespace and
        return ``None``. If ``False``, return a :class:`ScverseRFUResult`.
        DataFrame inputs always return :class:`RFUTableResult` and are not
        mutated.
    key_added
        Observation-aligned Awkward chain-result key.
    cell_summary
        Add explicit cell-level ``obs`` columns. Disabled by default because
        cells may contain multiple eligible TRB chains.
    summary_policy
        Explicit multi-chain selection policy. The default is ambiguity-aware
        and emits no cell-level RFU when eligible chains disagree.

    Returns
    -------
    RFUTableResult, ScverseRFUResult, or None
        Table output, explicit native output, or ``None`` for in-place native
        annotation.

    Notes
    -----
    RFU identity remains exact-CDR3 based. Optional Scirpy chain indices affect
    only the explicit ``primary_vdj`` cell-summary policy, never chain-level
    RFU assignments.
    """
    execution = {
        "rfu_dir": rfu_dir,
        "mode": mode,
        "threshold": threshold,
        "deduplicate": deduplicate,
        "chunk_size": chunk_size,
        "max_workers": max_workers,
        "executor": executor,
        "resume": resume,
        "force_recompute": force_recompute,
        "wrapper_r_path": wrapper_r_path,
        "rscript_bin": rscript_bin,
        "extra_r_args": extra_r_args,
        "workdir": workdir,
    }
    if isinstance(data, pd.DataFrame):
        return call_rfu_table(data, **execution)
    return assign_rfu_native(
        data,
        inplace=inplace,
        key_added=key_added,
        airr_mod=airr_mod,
        airr_key=airr_key,
        chain_idx_key=chain_idx_key,
        cell_summary=cell_summary,
        summary_policy=summary_policy,
        summary_key_added=summary_key_added,
        **execution,
    )


def call_rfu(
    adata: Any,
    *,
    backend: str = "rfu_repo",
    rfu_dir: PathLike | None = None,
    mode: str = "standard",
    threshold: float = 0.6,
    deduplicate: bool = True,
    chunk_size: int | None = None,
    max_workers: int = 1,
    executor: str = "process",
    resume: bool = True,
    force_recompute: bool = False,
    chain: str = "TRB",
    airr_key: str = "airr",
    prefer_productive: bool = True,
    wrapper_r_path: PathLike = "r/run_rfu_repo.R",
    rscript_bin: str = "Rscript",
    extra_r_args: Sequence[str] | None = None,
    workdir: PathLike | None = None,
    out_key: str = "rfu",
) -> pd.DataFrame:
    """
    Unified RFU calling entrypoint.

    Parameters
    ----------
    backend
        Currently supported:
          - "rfu_repo": call upstream RFU repo via r/run_rfu_repo.R
    rfu_dir
        Path to upstream RFU checkout. Falls back to RFU_DIR when omitted.
    mode
        ``"standard"`` (default) uses public ``AssignRFUs()``. ``"map_aware"``
        explicitly requests the optional ``AssignRFUs_with_map()`` capability.
    threshold
        RFU correlation threshold used for pass status and upstream ``N``.
    deduplicate
        Query each exact eligible CDR3 once, then restore row multiplicity.
    chunk_size
        Number of unique eligible CDR3 queries per chunk. ``None`` keeps
        the existing single-call behavior.
    max_workers
        Independent chunk workers. The default of one preserves serial execution.
    executor
        ``"process"`` or ``"thread"`` when chunk parallelism is enabled.
    resume
        Reuse only fully validated completed chunks when chunking is enabled.
    force_recompute
        Recompute every chunk even if a valid cache exists. Takes precedence
        over ``resume``.
    """
    backend = backend.lower().strip()
    if isinstance(extra_r_args, str):
        raise TypeError("extra_r_args must be a sequence of strings, not a single string.")

    if backend == "rfu_repo":
        return call_rfu_repo(
            adata,
            rfu_dir=rfu_dir,
            mode=mode,
            threshold=threshold,
            deduplicate=deduplicate,
            chunk_size=chunk_size,
            max_workers=max_workers,
            executor=executor,
            resume=resume,
            force_recompute=force_recompute,
            chain=chain,
            airr_key=airr_key,
            prefer_productive=prefer_productive,
            wrapper_r_path=wrapper_r_path,
            rscript_bin=rscript_bin,
            extra_r_args=extra_r_args,
            workdir=workdir,
            out_key=out_key,
        )

    raise ValueError(f"Unknown backend: {backend}. Supported: 'rfu_repo'")


__all__ = [
    "AntigenContextResult",
    "AntigenPermutationResult",
    "CohortHarmonizationResult",
    "ComparatorRepresentation",
    "FrozenRFUReference",
    "HeldOutValidationManifest",
    "LongitudinalCompartmentResult",
    "LongitudinalDesign",
    "LongitudinalDynamicsResult",
    "LongitudinalResamplingResult",
    "RFULongitudinalResult",
    "RFUOverlapResult",
    "RFUPseudobulkResult",
    "RFUTableResult",
    "ScverseRFUResult",
    "StabilityBenchmarkResult",
    "TransferCohortResult",
    "VDJdbEvidenceSummary",
    "VDJdbReference",
    "aggregate_rfu",
    "assign_rfu",
    "annotate_vdjdb",
    "benchmark_representation_stability",
    "bootstrap_longitudinal_statistic",
    "call_rfu",
    "call_rfu_repo",
    "call_rfu_table",
    "compare_antigen_groupings",
    "concat_scrfu",
    "create_heldout_validation_manifest",
    "deterministic_subsample",
    "donor_retrieval",
    "donor_leave_one_out",
    "global_antigen_coherence",
    "harmonize_cohort_metadata",
    "list_comparators",
    "load_vdjdb_reference",
    "longitudinal_compartment_comparison",
    "longitudinal_similarity",
    "multinomial_abundance_resample",
    "normalize_vdjdb_cdr3",
    "normalize_vdjdb_v_gene",
    "permute_longitudinal_labels",
    "reference_coverage",
    "register_comparator",
    "repertoire_metrics",
    "repertoire_representation",
    "rfu_antigen_abundance",
    "rfu_antigen_coherence",
    "rfu_antigen_permutation_test",
    "rfu_longitudinal_dynamics",
    "rfu_longitudinal_matrix",
    "rfu_metrics",
    "rfu_overlap",
    "rfu_phenotype_coupling",
    "rfu_pseudobulk",
    "rfu_sequence_matrix",
    "rfu_summary",
    "shuffle_input_order",
    "summarize_antigen_context",
    "summarize_longitudinal_similarity",
    "summarize_vdjdb_evidence",
    "threshold_sensitivity",
    "transfer_cohort",
    "validate_airr",
    "validate_completed_rfu_run",
    "validate_frozen_reference",
    "validate_longitudinal_design",
    "validate_scrfu_schema",
    "validate_vdjdb_reference",
]
