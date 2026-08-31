from __future__ import annotations

import hashlib
import json
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd

from ._version import __version__
from .pp import canonicalize_receptor_table, normalize_chain, normalize_productive, normalize_text
from .tl_rfu_repo import RFUTableResult, call_rfu_table

SCRFU_STORAGE_SCHEMA_VERSION = "1.0"
SCRFU_CHAIN_FIELDS = (
    "eligible",
    "rfu_id",
    "rfu_label",
    "rfu_score",
    "threshold_pass",
    "assignment_status",
    "unique_sequence_id",
    "locus",
    "chain_index",
    "reference_identity",
)

SummaryPolicy = Literal[
    "ambiguity_aware",
    "highest_score",
    "highest_threshold_score",
    "first_eligible",
    "primary_vdj",
]
PathLike = str | Path


@dataclass(frozen=True)
class ScverseRFUResult:
    """RFU assignment result for an AnnData or MuData AIRR modality.

    ``chain_records`` is an Awkward array with one record per source AIRR
    chain. It can be assigned directly to ``adata.obsm[key_added]``. The
    original canonical table outputs remain available through ``table_result``.
    """

    table_result: RFUTableResult
    chain_records: Any
    cell_summary: pd.DataFrame | None
    provenance: dict[str, Any]
    airr_mod: str
    airr_key: str
    key_added: str


def _require_awkward() -> Any:
    try:
        import awkward as ak
    except ImportError as exc:  # pragma: no cover - exercised in core-only installations
        raise ImportError(
            "Native Scirpy AIRR support requires Awkward Array. Install scRFU with "
            "the 'scirpy' extra: pip install '.[scirpy]'."
        ) from exc
    return ak


def _resolve_airr_modality(data: Any, airr_mod: str) -> tuple[Any, bool]:
    if hasattr(data, "mod"):
        if airr_mod not in data.mod:
            raise KeyError(f"MuData has no modality {airr_mod!r}.")
        return data.mod[airr_mod], True
    if hasattr(data, "obsm") and hasattr(data, "obs"):
        return data, False
    raise TypeError("Expected a pandas DataFrame, AnnData, or MuData object.")


def _airr_records(target: Any, airr_key: str) -> tuple[Any, list[list[dict[str, Any]]]]:
    if airr_key not in target.obsm:
        raise KeyError(f"AIRR modality has no obsm[{airr_key!r}] value.")
    ak = _require_awkward()
    airr = target.obsm[airr_key]
    if not isinstance(airr, ak.Array):
        raise TypeError(
            f"obsm[{airr_key!r}] must be the current Scirpy Awkward AIRR representation; "
            f"received {type(airr).__name__}. The compatibility DataFrame path remains "
            "available through scrfu.tl.call_rfu()."
        )
    if len(airr) != target.n_obs:
        raise ValueError("AIRR records are not aligned with the observation axis.")
    records = ak.to_list(airr)
    if any(not isinstance(chains, list) for chains in records):
        raise ValueError("AIRR data must contain one variable-length chain list per observation.")
    return airr, records


def _as_bool(value: Any) -> bool:
    normalized = normalize_productive(value)
    return bool(normalized) if not pd.isna(normalized) else False


def _extract_native_chains(
    target: Any,
    *,
    airr_key: str,
) -> tuple[Any, list[list[dict[str, Any]]], pd.DataFrame, dict[str, tuple[int, int]]]:
    airr, records = _airr_records(target, airr_key)
    obs_names = pd.Index(target.obs_names.astype(str))
    if obs_names.has_duplicates:
        raise ValueError("AnnData AIRR observations must have unique observation names.")

    rows: list[dict[str, Any]] = []
    positions: dict[str, tuple[int, int]] = {}
    for obs_position, (cell_id, chains) in enumerate(zip(obs_names, records, strict=True)):
        for chain_position, record in enumerate(chains):
            if not isinstance(record, dict):
                raise ValueError("Every AIRR chain must be a record/dictionary.")
            locus = normalize_chain(record.get("locus", record.get("chain")))
            productive = _as_bool(record.get("productive"))
            cdr3aa = normalize_text(record.get("junction_aa", record.get("cdr3aa")))
            is_trb = not pd.isna(locus) and locus == "TRB"
            if not is_trb or not productive or pd.isna(cdr3aa):
                continue
            input_row_id = f"airr_{obs_position:08d}_{chain_position:04d}"
            positions[input_row_id] = (obs_position, chain_position)
            rows.append(
                {
                    "input_row_id": input_row_id,
                    "cell_id": cell_id,
                    "chain": "TRB",
                    "cdr3aa": cdr3aa,
                    "v_call": record.get("v_call"),
                    "productive": True,
                    "source_adapter": "scirpy_airr_awkward",
                    "source_row_id": f"{cell_id}:{chain_position}",
                }
            )
    columns = [
        "input_row_id",
        "cell_id",
        "chain",
        "cdr3aa",
        "v_call",
        "productive",
        "source_adapter",
        "source_row_id",
    ]
    frame = canonicalize_receptor_table(pd.DataFrame(rows, columns=columns))
    return airr, records, frame, positions


def _empty_table_result() -> RFUTableResult:
    return RFUTableResult(
        per_sequence=pd.DataFrame(),
        per_row=pd.DataFrame(),
        mapping=pd.DataFrame(),
        provenance={
            "scrfu_version": __version__,
            "backend": "not_invoked_no_productive_trb",
            "selected_receptor_row_count": 0,
            "table_level_api": True,
        },
    )


def _portable_reference_identity(provenance: dict[str, Any]) -> str:
    identity = {
        key: provenance.get(key)
        for key in (
            "backend",
            "backend_mode",
            "rfu_r_sha256",
            "trimer_rdata_sha256",
            "km5000_rdata_sha256",
            "wrapper_r_sha256",
            "rfu_threshold",
        )
        if provenance.get(key) is not None
    }
    if not identity or provenance.get("backend") == "not_invoked_no_productive_trb":
        return "unconfigured:no-productive-trb"
    payload = json.dumps(identity, sort_keys=True, separators=(",", ":")).encode()
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def _python_scalar(value: Any) -> Any:
    if value is None or pd.isna(value):
        return None
    if isinstance(value, np.generic):
        return value.item()
    return value


def _typed_jagged(values: list[list[Any]], kind: str) -> Any:
    ak = _require_awkward()
    counts = np.asarray([len(row) for row in values], dtype=np.int64)
    flat = [item for row in values for item in row]
    valid = np.asarray([item is not None for item in flat], dtype=bool)
    if kind == "bool":
        content = np.asarray([False if item is None else bool(item) for item in flat], dtype=bool)
        array = ak.Array(content)
    elif kind == "int":
        content = np.asarray([0 if item is None else int(item) for item in flat], dtype=np.int64)
        array = ak.Array(content)
    elif kind == "float":
        content = np.asarray(
            [np.nan if item is None else float(item) for item in flat], dtype=float
        )
        array = ak.Array(content)
    elif kind == "str":
        content = ["" if item is None else str(item) for item in flat]
        array = ak.Array([""] + content)[1:]
    else:  # pragma: no cover - internal misuse guard
        raise ValueError(f"Unknown Awkward storage kind: {kind}")
    array = ak.mask(array, valid)
    return ak.unflatten(array, counts)


def _build_chain_records(
    source_records: list[list[dict[str, Any]]],
    table_result: RFUTableResult,
    positions: dict[str, tuple[int, int]],
    reference_identity: str,
) -> Any:
    by_position: dict[tuple[int, int], dict[str, Any]] = {}
    if not table_result.per_row.empty:
        if table_result.per_row["input_row_id"].duplicated().any():
            raise ValueError("RFU output contains duplicate input_row_id values.")
        returned = set(table_result.per_row["input_row_id"].astype(str))
        expected = set(positions)
        if returned != expected:
            raise ValueError(
                "RFU output identifiers do not exactly reconstruct native AIRR chain queries."
            )
        for row in table_result.per_row.to_dict(orient="records"):
            by_position[positions[str(row["input_row_id"])]] = row

    nested: dict[str, list[list[Any]]] = {field: [] for field in SCRFU_CHAIN_FIELDS}
    for obs_position, chains in enumerate(source_records):
        row_values = {field: [] for field in SCRFU_CHAIN_FIELDS}
        for chain_position, record in enumerate(chains):
            locus = normalize_chain(record.get("locus", record.get("chain")))
            productive = _as_bool(record.get("productive"))
            cdr3aa = normalize_text(record.get("junction_aa", record.get("cdr3aa")))
            assignment = by_position.get((obs_position, chain_position))
            if assignment is not None:
                eligible = assignment.get("eligibility_status") == "eligible"
                status = assignment.get("assignment_status", assignment.get("rfu_status"))
                values = {
                    "eligible": eligible,
                    "rfu_id": _python_scalar(assignment.get("rfu_id")),
                    "rfu_label": _python_scalar(assignment.get("rfu_label")),
                    "rfu_score": _python_scalar(assignment.get("rfu_score")),
                    "threshold_pass": _python_scalar(
                        assignment.get("rfu_pass_threshold", assignment.get("pass_thr"))
                    ),
                    "assignment_status": _python_scalar(status),
                    "unique_sequence_id": _python_scalar(assignment.get("unique_sequence_id")),
                }
            else:
                if pd.isna(locus) or locus != "TRB":
                    status = "non_target_locus"
                elif not productive:
                    status = "nonproductive_chain"
                elif pd.isna(cdr3aa):
                    status = "missing_sequence"
                else:  # pragma: no cover - guarded by exact output reconstruction
                    status = "not_queried"
                values = {
                    "eligible": False,
                    "rfu_id": None,
                    "rfu_label": None,
                    "rfu_score": None,
                    "threshold_pass": None,
                    "assignment_status": status,
                    "unique_sequence_id": None,
                }
            values.update(
                {
                    "locus": None if pd.isna(locus) else str(locus),
                    "chain_index": chain_position,
                    "reference_identity": reference_identity,
                }
            )
            for field in SCRFU_CHAIN_FIELDS:
                row_values[field].append(values[field])
        for field in SCRFU_CHAIN_FIELDS:
            nested[field].append(row_values[field])

    ak = _require_awkward()
    kinds = {
        "eligible": "bool",
        "rfu_id": "int",
        "rfu_label": "str",
        "rfu_score": "float",
        "threshold_pass": "bool",
        "assignment_status": "str",
        "unique_sequence_id": "str",
        "locus": "str",
        "chain_index": "int",
        "reference_identity": "str",
    }
    return ak.zip({field: _typed_jagged(nested[field], kinds[field]) for field in kinds})


def _primary_vdj_indices(target: Any, chain_idx_key: str) -> list[int | None]:
    ak = _require_awkward()
    if chain_idx_key not in target.obsm:
        raise KeyError(
            f"summary_policy='primary_vdj' requires obsm[{chain_idx_key!r}]. Run "
            "scirpy.pp.index_chains() first or choose another explicit summary policy."
        )
    indices = ak.to_list(target.obsm[chain_idx_key])
    result: list[int | None] = []
    for record in indices:
        vdj = record.get("VDJ", []) if isinstance(record, dict) else []
        result.append(next((int(value) for value in vdj if value is not None), None))
    return result


def _cell_summary(
    target: Any,
    chain_records: Any,
    *,
    policy: SummaryPolicy,
    chain_idx_key: str,
) -> pd.DataFrame:
    records = _require_awkward().to_list(chain_records)
    primary = _primary_vdj_indices(target, chain_idx_key) if policy == "primary_vdj" else None
    rows: list[dict[str, Any]] = []
    for obs_position, chains in enumerate(records):
        candidates = [
            chain for chain in chains if chain["eligible"] and chain["rfu_id"] is not None
        ]
        selected: dict[str, Any] | None = None
        summary_status = "no_eligible_assignment"
        if candidates:
            if policy == "ambiguity_aware":
                ids = {chain["rfu_id"] for chain in candidates}
                if len(ids) == 1:
                    selected = max(
                        candidates,
                        key=lambda item: (
                            float("-inf") if item["rfu_score"] is None else item["rfu_score"],
                            -item["chain_index"],
                        ),
                    )
                    summary_status = "unambiguous"
                else:
                    summary_status = "ambiguous_multiple_rfus"
            elif policy == "highest_score":
                selected = max(
                    candidates,
                    key=lambda item: (
                        float("-inf") if item["rfu_score"] is None else item["rfu_score"],
                        -item["chain_index"],
                    ),
                )
                summary_status = "selected_highest_score"
            elif policy == "highest_threshold_score":
                passed = [chain for chain in candidates if chain["threshold_pass"]]
                if passed:
                    selected = max(
                        passed,
                        key=lambda item: (
                            float("-inf") if item["rfu_score"] is None else item["rfu_score"],
                            -item["chain_index"],
                        ),
                    )
                    summary_status = "selected_highest_threshold_score"
                else:
                    summary_status = "no_threshold_qualified_assignment"
            elif policy == "first_eligible":
                selected = min(candidates, key=lambda item: item["chain_index"])
                summary_status = "selected_first_eligible"
            elif policy == "primary_vdj":
                index = primary[obs_position] if primary is not None else None
                selected = next(
                    (chain for chain in candidates if chain["chain_index"] == index), None
                )
                summary_status = (
                    "selected_primary_vdj" if selected is not None else "primary_vdj_unassigned"
                )
        rows.append(
            {
                "rfu_id": None if selected is None else selected["rfu_id"],
                "rfu_label": None if selected is None else selected["rfu_label"],
                "rfu_score": None if selected is None else selected["rfu_score"],
                "threshold_pass": None if selected is None else selected["threshold_pass"],
                "assignment_status": None if selected is None else selected["assignment_status"],
                "chain_index": None if selected is None else selected["chain_index"],
                "summary_status": summary_status,
            }
        )
    summary = pd.DataFrame(rows, index=target.obs_names.copy())
    summary["rfu_id"] = pd.array(summary["rfu_id"], dtype="Int64")
    summary["chain_index"] = pd.array(summary["chain_index"], dtype="Int64")
    summary["rfu_score"] = pd.to_numeric(summary["rfu_score"], errors="coerce")
    summary["threshold_pass"] = pd.array(summary["threshold_pass"], dtype="boolean")
    for column in ("rfu_label", "assignment_status", "summary_status"):
        summary[column] = summary[column].astype("string")
    return summary


def _native_provenance(
    table_result: RFUTableResult,
    *,
    airr_mod: str,
    airr_key: str,
    key_added: str,
    threshold: float,
    mode: str,
    summary_policy: str,
    reference_identity: str,
) -> dict[str, Any]:
    artifact_hashes = {
        key: str(table_result.provenance[key])
        for key in (
            "rfu_r_sha256",
            "trimer_rdata_sha256",
            "km5000_rdata_sha256",
            "wrapper_r_sha256",
        )
        if table_result.provenance.get(key)
    }
    return {
        "schema_version": SCRFU_STORAGE_SCHEMA_VERSION,
        "software_version": __version__,
        "reference_identity": reference_identity,
        "artifact_hashes": artifact_hashes,
        "assignment_mode": mode,
        "threshold": float(threshold),
        "target_locus": "TRB",
        "chain_selection": "all productive TRB chains with non-missing junction_aa",
        "summary_policy": summary_policy,
        "input_adapter": "scirpy_airr_awkward",
        "airr_mod": airr_mod,
        "airr_key": airr_key,
        "chain_result_key": key_added,
        "runtime_paths_stored": False,
    }


def assign_rfu_native(
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
    airr_mod: str = "airr",
    airr_key: str = "airr",
    chain_idx_key: str = "chain_indices",
    key_added: str = "scrfu",
    cell_summary: bool = False,
    summary_policy: SummaryPolicy = "ambiguity_aware",
    summary_key_added: str | None = None,
    wrapper_r_path: PathLike = "r/run_rfu_repo.R",
    rscript_bin: str = "Rscript",
    extra_r_args: Sequence[str] | None = None,
    workdir: PathLike | None = None,
    inplace: bool = True,
) -> ScverseRFUResult | None:
    """Assign RFUs to every eligible TRB chain in a Scirpy AIRR object.

    This function never reads or materializes ``X``. All source AIRR chains are
    retained in their original order, including non-TRB and nonproductive
    chains. RFU calculations delegate to the canonical table-level backend.
    """
    if summary_policy not in {
        "ambiguity_aware",
        "highest_score",
        "highest_threshold_score",
        "first_eligible",
        "primary_vdj",
    }:
        raise ValueError(f"Unknown cell-summary policy: {summary_policy!r}.")
    target, _ = _resolve_airr_modality(data, airr_mod)
    _, source_records, receptors, positions = _extract_native_chains(target, airr_key=airr_key)
    if receptors.empty:
        table_result = _empty_table_result()
    else:
        table_result = call_rfu_table(
            receptors,
            rfu_dir=rfu_dir,
            mode=mode,
            threshold=threshold,
            deduplicate=deduplicate,
            chunk_size=chunk_size,
            max_workers=max_workers,
            executor=executor,
            resume=resume,
            force_recompute=force_recompute,
            workdir=workdir,
            wrapper_r_path=wrapper_r_path,
            rscript_bin=rscript_bin,
            extra_r_args=extra_r_args,
        )
    reference_identity = _portable_reference_identity(table_result.provenance)
    chain_records = _build_chain_records(
        source_records, table_result, positions, reference_identity
    )
    summary = (
        _cell_summary(
            target,
            chain_records,
            policy=summary_policy,
            chain_idx_key=chain_idx_key,
        )
        if cell_summary
        else None
    )
    provenance = _native_provenance(
        table_result,
        airr_mod=airr_mod,
        airr_key=airr_key,
        key_added=key_added,
        threshold=threshold,
        mode=mode,
        summary_policy=summary_policy if cell_summary else "none",
        reference_identity=reference_identity,
    )
    result = ScverseRFUResult(
        table_result=table_result,
        chain_records=chain_records,
        cell_summary=summary,
        provenance=provenance,
        airr_mod=airr_mod,
        airr_key=airr_key,
        key_added=key_added,
    )
    if not inplace:
        return result

    target.obsm[key_added] = chain_records
    namespace = dict(target.uns.get("scrfu", {}))
    namespace.update(provenance)
    target.uns["scrfu"] = namespace
    if summary is not None:
        prefix = summary_key_added or key_added
        for column in summary:
            values = summary[column]
            # Object-backed text avoids requiring AnnData's opt-in nullable-string
            # writer and remains readable by older supported AnnData releases.
            if isinstance(values.dtype, pd.StringDtype):
                values = values.astype(object).where(values.notna(), None)
            target.obs[f"{prefix}_{column}"] = values
    validate_scrfu_schema(target, airr_key=airr_key, key_added=key_added)
    return None


def validate_scrfu_schema(
    data: Any,
    *,
    airr_mod: str = "airr",
    airr_key: str = "airr",
    key_added: str = "scrfu",
) -> dict[str, Any]:
    """Validate native scRFU chain alignment, fields, and portable provenance."""
    target, _ = _resolve_airr_modality(data, airr_mod)
    ak = _require_awkward()
    airr, _ = _airr_records(target, airr_key)
    if key_added not in target.obsm:
        raise KeyError(f"AIRR modality has no obsm[{key_added!r}] scRFU result.")
    result = target.obsm[key_added]
    if not isinstance(result, ak.Array):
        raise TypeError(f"obsm[{key_added!r}] is not an Awkward array.")
    if len(result) != len(airr):
        raise ValueError("scRFU results are not aligned with the AIRR observation axis.")
    airr_counts = ak.to_list(ak.num(airr, axis=1))
    result_counts = ak.to_list(ak.num(result, axis=1))
    if airr_counts != result_counts:
        raise ValueError("scRFU chain counts do not exactly match AIRR chain counts.")
    missing_fields = sorted(set(SCRFU_CHAIN_FIELDS).difference(ak.fields(result)))
    if missing_fields:
        raise ValueError(f"scRFU chain records are missing fields: {missing_fields}")
    namespace = target.uns.get("scrfu")
    if not isinstance(namespace, dict):
        raise ValueError("uns['scrfu'] portable provenance is missing.")
    if namespace.get("schema_version") != SCRFU_STORAGE_SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported scRFU storage schema version: {namespace.get('schema_version')!r}."
        )
    reference_identity = namespace.get("reference_identity")
    stored_references = {
        value
        for row in ak.to_list(result["reference_identity"])
        for value in row
        if value is not None
    }
    if stored_references and stored_references != {reference_identity}:
        raise ValueError("Chain records and provenance have incompatible RFU references.")
    return {
        "status": "ok",
        "schema_version": SCRFU_STORAGE_SCHEMA_VERSION,
        "observation_count": len(result),
        "chain_count": int(sum(result_counts)),
        "reference_identity": reference_identity,
        "airr_key": airr_key,
        "key_added": key_added,
    }


def concat_scrfu(
    objects: list[Any],
    *,
    key_added: str = "scrfu",
    **kwargs: Any,
) -> Any:
    """Concatenate compatible scRFU-annotated AnnData objects safely.

    All inputs must use the same storage schema and frozen RFU reference.
    This explicit guard prevents ordinary ``anndata.concat`` from silently
    dropping incompatible provenance through its ``uns_merge`` policy.
    """
    if not objects:
        raise ValueError("At least one AnnData object is required.")
    reports = [validate_scrfu_schema(obj, key_added=key_added) for obj in objects]
    identities = {report["reference_identity"] for report in reports}
    if len(identities) != 1:
        raise ValueError("Cannot concatenate scRFU objects with different RFU references.")
    versions = {obj.uns["scrfu"]["schema_version"] for obj in objects}
    if len(versions) != 1:
        raise ValueError("Cannot concatenate scRFU objects with different storage schemas.")
    import anndata as ad

    kwargs.setdefault("uns_merge", "same")
    combined = ad.concat(objects, **kwargs)
    validate_scrfu_schema(combined, key_added=key_added)
    return combined


__all__ = [
    "SCRFU_CHAIN_FIELDS",
    "SCRFU_STORAGE_SCHEMA_VERSION",
    "ScverseRFUResult",
    "assign_rfu_native",
    "concat_scrfu",
    "validate_scrfu_schema",
]
