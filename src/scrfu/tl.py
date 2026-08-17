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
    "StabilityBenchmarkResult",
    "TransferCohortResult",
    "VDJdbEvidenceSummary",
    "VDJdbReference",
    "aggregate_rfu",
    "annotate_vdjdb",
    "benchmark_representation_stability",
    "bootstrap_longitudinal_statistic",
    "call_rfu",
    "call_rfu_repo",
    "call_rfu_table",
    "compare_antigen_groupings",
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
    "validate_frozen_reference",
    "validate_longitudinal_design",
    "validate_vdjdb_reference",
]
