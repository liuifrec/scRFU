from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .attach import attach_rfu_results
from .backends.rfu_repo import RFURepoBackend, RFURepoPaths
from .extract import extract_trb_features
from .pp import canonicalize_receptor_table, validate_receptor_table

PathLike = str | Path


@dataclass(frozen=True)
class RFUTableResult:
    """Dataset-independent RFU outputs at sequence, receptor-row, and mapping levels."""

    per_sequence: pd.DataFrame
    per_row: pd.DataFrame
    mapping: pd.DataFrame
    provenance: dict[str, Any]


def call_rfu_table(
    receptors: pd.DataFrame,
    *,
    rfu_dir: PathLike | None = None,
    chain: str = "TRB",
    mode: str = "standard",
    threshold: float = 0.6,
    deduplicate: bool = True,
    chunk_size: int | None = None,
    max_workers: int = 1,
    executor: str = "process",
    resume: bool = True,
    force_recompute: bool = False,
    workdir: PathLike | None = None,
    wrapper_r_path: PathLike = "r/run_rfu_repo.R",
    rscript_bin: str = "Rscript",
    extra_r_args: Sequence[str] | None = None,
) -> RFUTableResult:
    """Call RFU directly from canonical receptor rows, without AnnData or dataset assumptions."""
    if isinstance(extra_r_args, str):
        raise TypeError("extra_r_args must be a sequence of strings, not a single string.")
    canonical = canonicalize_receptor_table(receptors)
    validate_receptor_table(canonical, strict=True)
    selected_chain = str(chain).strip().upper()
    if selected_chain != "TRB":
        raise ValueError(
            "The configured original RFU reference is TRB-specific; scRFU will not apply it "
            f"to chain {selected_chain!r}. A receptor-specific reference/backend is required."
        )
    selected = canonical.loc[canonical["chain"].eq(selected_chain)].copy()
    selected["trbv"] = selected["v_call"]

    backend = RFURepoBackend(
        rfu_dir=rfu_dir,
        mode=mode,
        wrapper_r_path=wrapper_r_path,
        rscript_bin=rscript_bin,
    )
    run = backend.run(
        selected,
        threshold=threshold,
        deduplicate=deduplicate,
        chunk_size=chunk_size,
        max_workers=max_workers,
        executor=executor,
        resume=resume,
        force_recompute=force_recompute,
        extra_args=extra_r_args,
        workdir=workdir,
    )
    per_row = run.df.copy()
    if "rfu_id_nearest" not in per_row:
        per_row["rfu_id_nearest"] = per_row["rfu_id"]
    if "rfu_label_nearest" not in per_row:
        per_row["rfu_label_nearest"] = per_row["rfu_label"]
    if "rfu_pass_threshold" not in per_row:
        per_row["rfu_pass_threshold"] = per_row["pass_thr"].astype("boolean")
    if "reference_coverage_status" not in per_row:
        eligible_mask = per_row["eligibility_status"].eq("eligible")
        assigned_mask = eligible_mask & per_row["rfu_id"].notna()
        passed_mask = per_row["pass_thr"].astype("boolean").fillna(False)
        coverage = pd.Series("ineligible_sequence", index=per_row.index, dtype="string")
        coverage.loc[per_row["cdr3aa"].isna()] = "missing_sequence"
        coverage.loc[eligible_mask & ~assigned_mask] = "upstream_unassigned"
        coverage.loc[assigned_mask & ~passed_mask] = "below_threshold"
        coverage.loc[assigned_mask & passed_mask] = "covered"
        per_row["reference_coverage_status"] = coverage
    if "assignment_status" not in per_row:
        per_row["assignment_status"] = per_row["reference_coverage_status"].map(
            {
                "covered": "nearest_threshold_qualified",
                "below_threshold": "nearest_below_threshold",
                "ineligible_sequence": "ineligible_sequence",
                "upstream_unassigned": "upstream_unassigned",
                "missing_sequence": "missing_sequence",
            }
        )
    mapping_columns = [
        column
        for column in (
            "input_row_id",
            "unique_sequence_id",
            "cell_id",
            "chain",
            "cdr3aa",
            "v_call",
            "source_adapter",
            "source_row_id",
            "eligibility_status",
            "reference_coverage_status",
            "assignment_status",
        )
        if column in per_row
    ]
    mapping = per_row.loc[:, mapping_columns].copy()
    eligible = per_row.loc[per_row["eligibility_status"].eq("eligible")].copy()
    if eligible.empty:
        per_sequence = pd.DataFrame(
            columns=[
                "unique_sequence_id",
                "cdr3aa",
                "query_v_call",
                "multiplicity",
                "rfu_id",
                "rfu_label",
                "rfu_score",
                "pass_thr",
                "rfu_status",
                "rfu_id_nearest",
                "rfu_label_nearest",
                "rfu_pass_threshold",
                "reference_coverage_status",
                "assignment_status",
            ]
        )
    else:
        summary = eligible.groupby("unique_sequence_id", sort=False, as_index=False).agg(
            cdr3aa=("cdr3aa", "first"),
            query_v_call=("v_call", "first"),
            multiplicity=("input_row_id", "size"),
        )
        assignments = eligible.drop_duplicates("unique_sequence_id", keep="first").loc[
            :,
            [
                "unique_sequence_id",
                "rfu_id",
                "rfu_label",
                "rfu_score",
                "pass_thr",
                "rfu_status",
                "rfu_id_nearest",
                "rfu_label_nearest",
                "rfu_pass_threshold",
                "reference_coverage_status",
                "assignment_status",
            ],
        ]
        per_sequence = summary.merge(
            assignments, on="unique_sequence_id", how="left", sort=False, validate="one_to_one"
        )
    provenance = {
        **backend.provenance_dict(),
        **run.metadata,
        "chain": selected_chain,
        "input_receptor_row_count": len(canonical),
        "selected_receptor_row_count": len(selected),
        "table_level_api": True,
    }
    return RFUTableResult(per_sequence, per_row, mapping, provenance)


def call_rfu_repo(
    adata: Any,
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
    Call upstream RFU (https://github.com/s175573/RFU) through scrfu's wrapper and attach results to AnnData.

    Writes:
      - adata.obs["trb_cdr3aa"], adata.obs["trbv"]
      - adata.obs["rfu_label"], adata.obs["rfu_score"]
      - adata.uns["scrfu"] provenance (merged)

    Standard mode is canonical and requires only public AssignRFUs(). The
    optional map-aware mode must be requested explicitly.

    Returns one result row per extracted feature row, including stable input and
    query identifiers, eligibility/assignment status, label, and score.
    """
    if isinstance(extra_r_args, str):
        raise TypeError("extra_r_args must be a sequence of strings, not a single string.")
    if str(chain).strip().upper() != "TRB":
        raise ValueError(
            "The configured original RFU reference is TRB-specific; another receptor chain "
            "requires a receptor-specific reference/backend."
        )

    # Preserve the historical configuration-error precedence before input extraction.
    RFURepoPaths.resolve(rfu_dir)

    features = extract_trb_features(
        adata,
        airr_key=airr_key,
        chain=chain,
        prefer_productive=prefer_productive,
    )

    canonical = pd.DataFrame(
        {
            "input_row_id": range(len(features)),
            "cell_id": features["cell_id"],
            "chain": str(chain).upper(),
            "cdr3aa": features["cdr3aa"],
            "v_call": features["trbv"],
            "productive": True,
            "source_adapter": "anndata_compat",
            "source_row_id": features.index.astype(str),
        }
    )
    canonical = canonicalize_receptor_table(canonical, source_adapter="anndata_compat")
    table_run = call_rfu_table(
        canonical,
        rfu_dir=rfu_dir,
        chain=chain,
        mode=mode,
        threshold=threshold,
        deduplicate=deduplicate,
        chunk_size=chunk_size,
        max_workers=max_workers,
        executor=executor,
        resume=resume,
        force_recompute=force_recompute,
        extra_r_args=extra_r_args,
        workdir=workdir,
        wrapper_r_path=wrapper_r_path,
        rscript_bin=rscript_bin,
    )

    # Attach results + provenance
    attach_rfu_results(
        adata,
        features=features,
        rfu_df=table_run.per_row,
        provenance={
            **table_run.provenance,
            "airr_key": airr_key,
            "out_key": out_key,
        },
        out_key=out_key,
    )
    return table_run.per_row


__all__ = ["RFUTableResult", "call_rfu_repo", "call_rfu_table"]
