from __future__ import annotations

import hashlib
import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ._version import __version__
from .pp import normalize_chain, normalize_text
from .repertoire import analysis_frame, shannon

VDJDB_REFERENCE_SCHEMA_VERSION = "1.0"
TCR_CHAINS = frozenset({"TRA", "TRB", "TRG", "TRD"})
AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")

REFERENCE_COLUMNS = (
    "reference_row_id",
    "cdr3aa",
    "v_call",
    "j_call",
    "chain",
    "paired_receptor_id",
    "paired_cdr3aa",
    "paired_v_call",
    "epitope",
    "antigen_gene",
    "antigen_species",
    "mhc",
    "mhc_class",
    "evidence_score",
    "publication_id",
    "database_row_id",
)

_ALIASES: dict[str, tuple[str, ...]] = {
    "cdr3aa": ("cdr3aa", "cdr3", "cdr3_aa", "cdr3 amino acid"),
    "v_call": ("v_call", "v", "v gene", "v_gene", "v.segm", "v segm"),
    "j_call": ("j_call", "j", "j gene", "j_gene", "j.segm", "j segm"),
    "chain": ("chain", "gene", "locus", "receptor chain"),
    "paired_receptor_id": ("paired_receptor_id", "complex.id", "complex_id", "pair_id"),
    "paired_cdr3aa": ("paired_cdr3aa", "paired_cdr3", "cdr3_pair"),
    "paired_v_call": ("paired_v_call", "paired_v", "v_pair"),
    "epitope": ("epitope", "epitope sequence", "antigen epitope"),
    "antigen_gene": ("antigen_gene", "epitope gene", "antigen gene"),
    "antigen_species": ("antigen_species", "epitope species", "species"),
    "mhc": ("mhc", "mhc allele", "mhc a", "hla"),
    "mhc_class": ("mhc_class", "mhc class"),
    "evidence_score": ("evidence_score", "score", "vdjdb score"),
    "publication_id": (
        "publication_id",
        "reference",
        "reference.id",
        "reference id",
        "publication",
        "pmid",
    ),
    "database_row_id": ("database_row_id", "row_id", "vdjdb_id", "id"),
}


@dataclass(frozen=True)
class VDJdbReference:
    table: pd.DataFrame
    provenance: dict[str, Any]
    validation: dict[str, Any]


@dataclass(frozen=True)
class VDJdbEvidenceSummary:
    sequence_summary: pd.DataFrame
    row_summary: pd.DataFrame


@dataclass(frozen=True)
class AntigenPermutationResult:
    observed: float
    permutation_values: np.ndarray
    null_mean: float
    null_std: float
    empirical_upper_tail_probability: float
    z_score: float
    parameters: dict[str, Any]


@dataclass(frozen=True)
class AntigenContextResult:
    sequence_prevalence: pd.DataFrame
    rfu_antigen_recurrence: pd.DataFrame
    phenotype_abundance: pd.DataFrame
    phenotype_evidence_tiers: pd.DataFrame
    phenotype_ambiguity: pd.DataFrame
    join_qc: dict[str, Any]


def _header(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(value).strip().lower()).strip()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_vdjdb_cdr3(value: Any, *, validate: bool = True) -> Any:
    """Trim and uppercase a CDR3 without altering amino-acid characters."""
    normalized = normalize_text(value)
    if pd.isna(normalized):
        return pd.NA
    result = str(normalized).upper()
    invalid = sorted(set(result).difference(AMINO_ACIDS))
    if validate and invalid:
        raise ValueError(f"CDR3 contains invalid amino-acid characters: {invalid}")
    return result


def normalize_vdjdb_v_gene(value: Any, *, mode: str = "strip_allele") -> Any:
    """Normalize V genes exactly or by removing only the terminal allele suffix."""
    if mode not in {"exact", "strip_allele"}:
        raise ValueError("v_gene_mode must be 'exact' or 'strip_allele'.")
    normalized = normalize_text(value)
    if pd.isna(normalized):
        return pd.NA
    result = str(normalized).upper()
    return re.sub(r"\*[A-Z0-9]+$", "", result) if mode == "strip_allele" else result


def _resolve_columns(
    frame: pd.DataFrame, overrides: Mapping[str, str] | None
) -> dict[str, Any | None]:
    normalized = {_header(column): column for column in frame.columns}
    mapping: dict[str, Any | None] = {}
    overrides = dict(overrides or {})
    unknown = sorted(set(overrides).difference(_ALIASES))
    if unknown:
        raise ValueError(f"Unknown canonical VDJdb column overrides: {unknown}")
    for canonical, aliases in _ALIASES.items():
        if canonical in overrides:
            source = overrides[canonical]
            if source not in frame:
                raise ValueError(f"VDJdb column override {canonical!r} -> {source!r} is missing.")
            mapping[canonical] = source
        else:
            mapping[canonical] = next(
                (normalized[_header(alias)] for alias in aliases if _header(alias) in normalized),
                None,
            )
    return mapping


def validate_vdjdb_reference(
    reference: VDJdbReference | pd.DataFrame, *, strict: bool = True
) -> dict[str, Any]:
    """Validate the canonical local-reference table without collapsing annotations."""
    table = reference.table if isinstance(reference, VDJdbReference) else reference
    if not isinstance(table, pd.DataFrame):
        raise TypeError("VDJdb reference must be a VDJdbReference or pandas DataFrame.")
    errors: list[str] = []
    if "cdr3aa" not in table:
        errors.append("missing required CDR3 field 'cdr3aa'")
        missing = len(table)
        invalid_rows: list[int] = []
    else:
        missing = int(table["cdr3aa"].isna().sum())
        invalid_rows = [
            int(index)
            for index, value in table["cdr3aa"].items()
            if not pd.isna(value) and set(str(value)).difference(AMINO_ACIDS)
        ]
        if missing:
            errors.append(f"cdr3aa contains {missing} missing values")
        if invalid_rows:
            errors.append(
                f"cdr3aa contains invalid amino-acid characters at rows {invalid_rows[:10]}"
            )
    if "reference_row_id" not in table:
        errors.append("missing reference_row_id")
        duplicate_ids: list[str] = []
    else:
        duplicate_ids = (
            table.loc[table["reference_row_id"].duplicated(keep=False), "reference_row_id"]
            .astype(str)
            .drop_duplicates()
            .tolist()
        )
        if duplicate_ids:
            errors.append(f"duplicate reference_row_id values: {duplicate_ids[:10]}")
    annotation_columns = [
        column
        for column in ("epitope", "antigen_gene", "antigen_species", "mhc", "mhc_class")
        if column in table and table[column].notna().any()
    ]
    duplicate_database_rows = int(
        table.drop(columns=["reference_row_id"], errors="ignore").duplicated(keep=False).sum()
    )
    report = {
        "schema_version": VDJDB_REFERENCE_SCHEMA_VERSION,
        "status": "ok" if not errors else "error",
        "row_count": len(table),
        "unique_cdr3_count": int(table["cdr3aa"].nunique()) if "cdr3aa" in table else 0,
        "missing_cdr3_count": missing,
        "invalid_cdr3_rows": invalid_rows,
        "duplicate_reference_row_ids": duplicate_ids,
        "duplicated_database_row_count": duplicate_database_rows,
        "available_chains": sorted(
            table.get("chain", pd.Series(dtype="string")).dropna().astype(str).unique()
        ),
        "available_antigen_fields": annotation_columns,
        "errors": errors,
    }
    if strict and errors:
        raise ValueError("Invalid VDJdb reference: " + "; ".join(errors))
    return report


def load_vdjdb_reference(
    source: str | Path | pd.DataFrame,
    *,
    release_label: str,
    source_url: str | None = None,
    expected_sha256: str | None = None,
    column_mappings: Mapping[str, str] | None = None,
    sep: str | None = None,
    strict: bool = True,
) -> VDJdbReference:
    """Load a user-supplied, explicitly version-labelled local VDJdb table."""
    if not str(release_label).strip():
        raise ValueError("An explicit non-empty VDJdb release_label is required.")
    if isinstance(source, pd.DataFrame):
        raw = source.copy()
        payload = raw.to_csv(index=False, lineterminator="\n").encode()
        digest = hashlib.sha256(payload).hexdigest()
        source_path = None
        source_name = "<dataframe>"
        size = len(payload)
    else:
        path = Path(source).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"VDJdb reference not found: {path}")
        digest = _sha256(path)
        source_path = str(path)
        source_name = path.name
        size = path.stat().st_size
        inferred = (
            "\t" if ".tsv" in path.name.lower() else "," if ".csv" in path.name.lower() else None
        )
        read_sep = sep if sep is not None else inferred
        raw = pd.read_csv(path, sep=read_sep, engine="python" if read_sep is None else "c")
    if expected_sha256 is not None and digest.lower() != expected_sha256.strip().lower():
        raise ValueError(f"VDJdb SHA256 mismatch: expected {expected_sha256}, observed {digest}.")
    mapping = _resolve_columns(raw, column_mappings)
    if mapping["cdr3aa"] is None:
        raise ValueError(
            f"VDJdb reference is missing a required CDR3 field. Columns: {list(raw.columns)}"
        )
    table = pd.DataFrame(index=raw.index)
    for canonical in _ALIASES:
        source_column = mapping[canonical]
        table[canonical] = raw[source_column] if source_column is not None else pd.NA
    table.insert(0, "reference_row_id", [f"vdjdb_{index:09d}" for index in range(len(table))])
    for column in ("cdr3aa", "paired_cdr3aa"):
        table[column] = (
            table[column]
            .map(lambda value: normalize_vdjdb_cdr3(value, validate=False))
            .astype("string")
        )
    for column in ("v_call", "paired_v_call", "j_call"):
        table[column] = table[column].map(normalize_text).astype("string")
    table["chain"] = table["chain"].map(normalize_chain).astype("string")
    for column in (
        "paired_receptor_id",
        "epitope",
        "antigen_gene",
        "antigen_species",
        "mhc",
        "mhc_class",
        "publication_id",
        "database_row_id",
    ):
        table[column] = table[column].map(normalize_text).astype("string")
    table["evidence_score"] = pd.to_numeric(table["evidence_score"], errors="coerce")
    table = table.loc[:, REFERENCE_COLUMNS]
    validation = validate_vdjdb_reference(table, strict=strict)
    provenance = {
        "reference_schema_version": VDJDB_REFERENCE_SCHEMA_VERSION,
        "release_label": str(release_label).strip(),
        "source_url": source_url,
        "source_filename": source_name,
        "source_path": source_path,
        "file_size": size,
        "sha256": digest,
        "load_timestamp": datetime.now(timezone.utc).isoformat(),
        "row_count": len(table),
        "unique_cdr3_count": int(table["cdr3aa"].nunique()),
        "available_receptor_chains": validation["available_chains"],
        "available_antigen_fields": validation["available_antigen_fields"],
        "column_mappings": {
            canonical: str(source_column) if source_column is not None else None
            for canonical, source_column in mapping.items()
        },
        "scrfu_version": __version__,
    }
    return VDJdbReference(table.reset_index(drop=True), provenance, validation)


def _prepare_queries(
    data: Any,
    *,
    chain: str | None,
    v_gene_mode: str,
    cdr3_col: str = "cdr3aa",
    v_col: str = "v_call",
    chain_col: str = "chain",
) -> pd.DataFrame:
    frame = analysis_frame(data).reset_index(drop=True)
    if cdr3_col not in frame:
        raise ValueError(f"VDJdb matching input is missing CDR3 column {cdr3_col!r}.")
    if "input_row_id" not in frame:
        frame["input_row_id"] = [f"row_{index:09d}" for index in range(len(frame))]
    if frame["input_row_id"].isna().any() or frame["input_row_id"].duplicated().any():
        raise ValueError("VDJdb matching requires unique, non-missing input_row_id values.")
    frame["_source_order"] = range(len(frame))
    frame["_cdr3_norm"] = frame[cdr3_col].map(normalize_vdjdb_cdr3).astype("string")
    frame["_v_norm"] = (
        frame[v_col]
        .map(lambda value: normalize_vdjdb_v_gene(value, mode=v_gene_mode))
        .astype("string")
        if v_col in frame
        else pd.Series(pd.NA, index=frame.index, dtype="string")
    )
    frame["_chain_norm"] = (
        frame[chain_col].map(normalize_chain).astype("string")
        if chain_col in frame
        else pd.Series(pd.NA, index=frame.index, dtype="string")
    )
    if chain is not None:
        selected_chain = normalize_chain(chain)
        frame = frame.loc[frame["_chain_norm"].eq(selected_chain).fillna(False)].copy()
    if "unique_sequence_id" not in frame:
        keys = list(zip(frame["_cdr3_norm"], frame["_v_norm"], frame["_chain_norm"], strict=True))
        lookup: dict[tuple[Any, ...], str] = {}
        ids: list[str] = []
        for key in keys:
            if key not in lookup:
                lookup[key] = f"query_{len(lookup):09d}"
            ids.append(lookup[key])
        frame["unique_sequence_id"] = ids
    consistency = frame.groupby("unique_sequence_id", observed=True)[
        ["_cdr3_norm", "_v_norm", "_chain_norm"]
    ].nunique(dropna=False)
    if (consistency > 1).any().any():
        raise ValueError("A unique_sequence_id maps to conflicting receptor features.")
    return frame


def annotate_vdjdb(
    receptors_or_rfu_results: Any,
    reference: VDJdbReference,
    *,
    match_mode: str = "cdr3_v",
    chain: str | None = "TRB",
    v_gene_mode: str = "strip_allele",
) -> pd.DataFrame:
    """Return cell/input-row-expanded exact evidence without mutating query results."""
    if not isinstance(reference, VDJdbReference):
        raise TypeError("reference must be loaded with load_vdjdb_reference().")
    if match_mode not in {"cdr3", "cdr3_v", "paired_exact"}:
        raise ValueError("match_mode must be 'cdr3', 'cdr3_v', or 'paired_exact'.")
    queries = _prepare_queries(receptors_or_rfu_results, chain=chain, v_gene_mode=v_gene_mode)
    refs = reference.table.copy()
    refs["_cdr3_norm"] = (
        refs["cdr3aa"]
        .map(lambda value: normalize_vdjdb_cdr3(value, validate=False))
        .astype("string")
    )
    refs["_v_norm"] = (
        refs["v_call"]
        .map(lambda value: normalize_vdjdb_v_gene(value, mode=v_gene_mode))
        .astype("string")
    )
    refs["_chain_norm"] = refs["chain"].map(normalize_chain).astype("string")
    if chain is not None:
        selected_chain = normalize_chain(chain)
        refs = refs.loc[refs["_chain_norm"].isna() | refs["_chain_norm"].eq(selected_chain)].copy()
    unique = queries.drop_duplicates("unique_sequence_id", keep="first").copy()
    keys = ["_cdr3_norm"]
    if match_mode in {"cdr3_v", "paired_exact"}:
        if unique["_v_norm"].isna().any():
            unique = unique.loc[unique["_v_norm"].notna()].copy()
        refs = refs.loc[refs["_v_norm"].notna()].copy()
        keys.append("_v_norm")
    if match_mode == "paired_exact":
        required = ["paired_cdr3aa"]
        missing = [
            column for column in required if column not in unique or unique[column].isna().all()
        ]
        if missing or refs["paired_cdr3aa"].isna().all():
            raise ValueError(
                "paired_exact requires genuine paired_cdr3aa fields in query and reference."
            )
        queries["_paired_cdr3_norm"] = (
            queries["paired_cdr3aa"].map(normalize_vdjdb_cdr3).astype("string")
        )
        paired_consistency = queries.groupby("unique_sequence_id", observed=True)[
            "_paired_cdr3_norm"
        ].nunique(dropna=False)
        if (paired_consistency > 1).any():
            raise ValueError("A unique_sequence_id maps to conflicting paired CDR3 sequences.")
        unique["_paired_cdr3_norm"] = (
            unique["paired_cdr3aa"].map(normalize_vdjdb_cdr3).astype("string")
        )
        refs["_paired_cdr3_norm"] = refs["paired_cdr3aa"].map(normalize_vdjdb_cdr3).astype("string")
        keys.append("_paired_cdr3_norm")
        if (
            "paired_v_call" in unique
            and unique["paired_v_call"].notna().any()
            and refs["paired_v_call"].notna().any()
        ):
            unique["_paired_v_norm"] = (
                unique["paired_v_call"]
                .map(lambda value: normalize_vdjdb_v_gene(value, mode=v_gene_mode))
                .astype("string")
            )
            refs["_paired_v_norm"] = (
                refs["paired_v_call"]
                .map(lambda value: normalize_vdjdb_v_gene(value, mode=v_gene_mode))
                .astype("string")
            )
            keys.append("_paired_v_norm")
    refs = refs.rename(
        columns={
            "cdr3aa": "matched_cdr3aa",
            "v_call": "matched_v_call",
            "chain": "matched_chain",
            "j_call": "matched_j_call",
            "paired_receptor_id": "reference_paired_receptor_id",
            "paired_cdr3aa": "matched_paired_cdr3aa",
            "paired_v_call": "matched_paired_v_call",
            "_chain_norm": "_reference_chain_norm",
        }
    )
    reference_payload = [
        "reference_row_id",
        "matched_cdr3aa",
        "matched_v_call",
        "matched_j_call",
        "matched_chain",
        "reference_paired_receptor_id",
        "matched_paired_cdr3aa",
        "matched_paired_v_call",
        "epitope",
        "antigen_gene",
        "antigen_species",
        "mhc",
        "mhc_class",
        "evidence_score",
        "publication_id",
        "database_row_id",
        "_reference_chain_norm",
        *keys,
    ]
    matched = unique.merge(
        refs[reference_payload], on=keys, how="inner", sort=False, validate="many_to_many"
    )
    matched = matched.loc[
        matched["_reference_chain_norm"].isna()
        | matched["_chain_norm"].isna()
        | matched["_reference_chain_norm"].eq(matched["_chain_norm"])
    ].copy()
    if matched.duplicated(["unique_sequence_id", "reference_row_id"]).any():
        raise ValueError("VDJdb matching produced duplicate sequence/reference pairs.")
    matched = matched.rename(
        columns={
            "_cdr3_norm": "query_cdr3aa",
            "_v_norm": "query_v_call",
            "_chain_norm": "query_chain",
            "paired_cdr3aa": "query_paired_cdr3aa",
            "paired_v_call": "query_paired_v_call",
        }
    )
    tier = {"cdr3": "cdr3_exact", "cdr3_v": "cdr3_v_exact", "paired_exact": "paired_exact"}[
        match_mode
    ]
    matched["match_mode"] = match_mode
    matched["evidence_tier"] = tier
    matched["reference_release"] = reference.provenance["release_label"]
    matched["reference_sha256"] = reference.provenance["sha256"]
    sequence_evidence_columns = [
        "unique_sequence_id",
        "query_cdr3aa",
        "query_v_call",
        "query_chain",
        "reference_row_id",
        "matched_cdr3aa",
        "matched_v_call",
        "matched_j_call",
        "matched_chain",
        "query_paired_cdr3aa",
        "query_paired_v_call",
        "reference_paired_receptor_id",
        "matched_paired_cdr3aa",
        "matched_paired_v_call",
        "epitope",
        "antigen_gene",
        "antigen_species",
        "mhc",
        "mhc_class",
        "evidence_score",
        "publication_id",
        "database_row_id",
        "match_mode",
        "evidence_tier",
        "reference_release",
        "reference_sha256",
    ]
    matched = matched.reindex(columns=sequence_evidence_columns)
    row_columns = [
        "unique_sequence_id",
        "input_row_id",
        "_source_order",
        *[
            column
            for column in ("cell_id", "source_row_id", "rfu_label", "pass_thr")
            if column in queries
        ],
    ]
    expected_rows = int(
        matched.groupby("unique_sequence_id", observed=True)
        .size()
        .mul(queries.groupby("unique_sequence_id", observed=True).size(), fill_value=0)
        .sum()
    )
    evidence = matched.merge(
        queries[row_columns],
        on="unique_sequence_id",
        how="left",
        sort=False,
        validate="many_to_many",
    )
    if (
        len(evidence) != expected_rows
        or evidence.duplicated(["input_row_id", "reference_row_id"]).any()
    ):
        raise RuntimeError("VDJdb evidence reconstruction produced an unexpected row expansion.")
    evidence = evidence.sort_values(["_source_order", "reference_row_id"], kind="stable")
    output = [
        "unique_sequence_id",
        "input_row_id",
        *[
            column
            for column in ("cell_id", "source_row_id", "rfu_label", "pass_thr")
            if column in evidence
        ],
        "query_cdr3aa",
        "query_v_call",
        "query_chain",
        "reference_row_id",
        "matched_cdr3aa",
        "matched_v_call",
        "matched_j_call",
        "matched_chain",
        "query_paired_cdr3aa",
        "query_paired_v_call",
        "reference_paired_receptor_id",
        "matched_paired_cdr3aa",
        "matched_paired_v_call",
        "epitope",
        "antigen_gene",
        "antigen_species",
        "mhc",
        "mhc_class",
        "evidence_score",
        "publication_id",
        "database_row_id",
        "match_mode",
        "evidence_tier",
        "reference_release",
        "reference_sha256",
    ]
    return evidence.reindex(columns=output).reset_index(drop=True)


def summarize_vdjdb_evidence(queries: Any, evidence: pd.DataFrame) -> VDJdbEvidenceSummary:
    """Summarize authoritative long evidence without selecting an ambiguous epitope."""
    rows = _prepare_queries(queries, chain=None, v_gene_mode="strip_allele")
    sequences = rows.drop_duplicates("unique_sequence_id", keep="first")
    pairs = evidence.drop_duplicates(["unique_sequence_id", "reference_row_id"]).copy()
    summaries: list[dict[str, Any]] = []
    for sequence_id in sequences["unique_sequence_id"]:
        subset = pairs.loc[pairs["unique_sequence_id"].eq(sequence_id)]
        epitopes = (
            sorted(subset["epitope"].dropna().astype(str).unique()) if "epitope" in subset else []
        )
        tiers = (
            sorted(subset["evidence_tier"].dropna().astype(str).unique())
            if "evidence_tier" in subset
            else []
        )
        score = pd.to_numeric(subset.get("evidence_score", pd.Series(dtype=float)), errors="coerce")
        summaries.append(
            {
                "unique_sequence_id": sequence_id,
                "has_vdjdb_evidence": not subset.empty,
                "evidence_record_count": len(subset),
                "distinct_epitope_count": len(epitopes),
                "dominant_epitope": epitopes[0] if len(epitopes) == 1 else pd.NA,
                "antigen_ambiguity": len(epitopes) > 1,
                "maximum_evidence_score": float(score.max())
                if score.notna().any()
                else float("nan"),
                "evidence_tiers": "|".join(tiers),
            }
        )
    sequence_summary = (
        sequences[["unique_sequence_id", "_cdr3_norm", "_v_norm", "_chain_norm"]]
        .rename(columns={"_cdr3_norm": "cdr3aa", "_v_norm": "v_call", "_chain_norm": "chain"})
        .merge(pd.DataFrame(summaries), on="unique_sequence_id", how="left", validate="one_to_one")
    )
    row_summary = rows.drop(columns=[column for column in rows if column.startswith("_")]).merge(
        sequence_summary.drop(columns=["cdr3aa", "v_call", "chain"]),
        on="unique_sequence_id",
        how="left",
        sort=False,
        validate="many_to_one",
    )
    if len(row_summary) != len(rows):
        raise RuntimeError("VDJdb row-summary reconstruction changed the input row count.")
    return VDJdbEvidenceSummary(
        sequence_summary.reset_index(drop=True), row_summary.reset_index(drop=True)
    )


def _assigned_unique(
    rfu_results: Any, *, assignment_policy: str, rfu_col: str = "rfu_label"
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if assignment_policy not in {"nearest", "threshold_pass"}:
        raise ValueError("assignment_policy must be 'nearest' or 'threshold_pass'.")
    rows = _prepare_queries(rfu_results, chain=None, v_gene_mode="strip_allele")
    if rfu_col not in rows:
        raise ValueError(f"RFU antigen analysis requires column {rfu_col!r}.")
    mask = rows[rfu_col].notna()
    if assignment_policy == "threshold_pass":
        if "pass_thr" not in rows:
            raise ValueError("threshold_pass policy requires pass_thr.")
        mask &= rows["pass_thr"].astype("boolean").fillna(False)
    rows = rows.loc[mask].copy()
    consistency = rows.groupby("unique_sequence_id", observed=True)[rfu_col].nunique()
    if (consistency > 1).any():
        raise ValueError("One unique sequence is assigned to multiple RFUs.")
    unique = rows.drop_duplicates("unique_sequence_id", keep="first")
    return rows, unique


def _sequence_antigens(
    evidence: pd.DataFrame, *, antigen_key: str, ambiguity_policy: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if ambiguity_policy not in {"fractional", "exclude_ambiguous", "multi_label"}:
        raise ValueError(
            "ambiguity_policy must be 'fractional', 'exclude_ambiguous', or 'multi_label'."
        )
    if antigen_key not in evidence:
        raise ValueError(f"Evidence is missing antigen field {antigen_key!r}.")
    pairs = evidence.drop_duplicates(["unique_sequence_id", "reference_row_id"])
    labels = pairs.loc[
        pairs[antigen_key].notna(), ["unique_sequence_id", antigen_key]
    ].drop_duplicates()
    counts = (
        labels.groupby("unique_sequence_id", observed=True)[antigen_key]
        .nunique()
        .rename("label_count")
    )
    labels = labels.merge(counts, on="unique_sequence_id", validate="many_to_one")
    labels["ambiguous"] = labels["label_count"].gt(1)
    if ambiguity_policy == "exclude_ambiguous":
        labels = labels.loc[~labels["ambiguous"]].copy()
    labels["label_weight"] = (
        1.0 / labels["label_count"] if ambiguity_policy == "fractional" else 1.0
    )
    sequence_status = counts.reset_index().assign(ambiguous=lambda x: x["label_count"].gt(1))
    return labels, sequence_status


def rfu_antigen_abundance(
    rfu_results: Any,
    evidence: pd.DataFrame,
    *,
    assignment_policy: str = "nearest",
    weighting: str = "unique_sequence",
    antigen_key: str = "epitope",
    ambiguity_policy: str = "fractional",
) -> pd.DataFrame:
    """Return long RFU-by-antigen abundance under an explicit ambiguity policy."""
    if weighting not in {"unique_sequence", "cell"}:
        raise ValueError("weighting must be 'unique_sequence' or 'cell'.")
    rows, unique = _assigned_unique(rfu_results, assignment_policy=assignment_policy)
    labels, _ = _sequence_antigens(
        evidence, antigen_key=antigen_key, ambiguity_policy=ambiguity_policy
    )
    if weighting == "cell":
        cell_identifier = "cell_id" if "cell_id" in rows else "input_row_id"
        units = (
            rows.groupby(["unique_sequence_id", "rfu_label"], observed=True)[cell_identifier]
            .nunique()
            .rename("unit_weight")
            .reset_index()
        )
    else:
        units = unique[["unique_sequence_id", "rfu_label"]].assign(unit_weight=1.0)
    joined = units.merge(labels, on="unique_sequence_id", how="inner", validate="one_to_many")
    joined["antigen_abundance"] = joined["unit_weight"] * joined["label_weight"]
    result = (
        joined.groupby(["rfu_label", antigen_key], observed=True)["antigen_abundance"]
        .sum()
        .reset_index()
    )
    totals = result.groupby("rfu_label", observed=True)["antigen_abundance"].transform("sum")
    result["within_rfu_proportion"] = result["antigen_abundance"] / totals
    result["weighting"] = weighting
    result["assignment_policy"] = assignment_policy
    result["ambiguity_policy"] = ambiguity_policy
    return result.sort_values(
        ["rfu_label", "antigen_abundance", antigen_key],
        ascending=[True, False, True],
        kind="stable",
    ).reset_index(drop=True)


def rfu_antigen_coherence(
    rfu_results: Any,
    evidence: pd.DataFrame,
    *,
    assignment_policy: str = "nearest",
    weighting: str = "unique_sequence",
    antigen_key: str = "epitope",
    min_matched_sequences: int = 2,
    ambiguity_policy: str = "fractional",
    sample_key: str | None = None,
    donor_key: str | None = None,
) -> pd.DataFrame:
    """Calculate descriptive RFU antigen-evidence coherence, never specificity."""
    if min_matched_sequences < 1:
        raise ValueError("min_matched_sequences must be positive.")
    rows, unique = _assigned_unique(rfu_results, assignment_policy=assignment_policy)
    pairs = evidence.drop_duplicates(["unique_sequence_id", "reference_row_id"])
    labels, status = _sequence_antigens(
        evidence, antigen_key=antigen_key, ambiguity_policy=ambiguity_policy
    )
    abundance = rfu_antigen_abundance(
        rfu_results,
        evidence,
        assignment_policy=assignment_policy,
        weighting=weighting,
        antigen_key=antigen_key,
        ambiguity_policy=ambiguity_policy,
    )
    labelled_ids = set(status["unique_sequence_id"])
    contributing_ids = set(labels["unique_sequence_id"])
    evidence_ids = set(pairs["unique_sequence_id"])
    output: list[dict[str, Any]] = []
    for rfu, rfu_unique in unique.groupby("rfu_label", sort=True, observed=True):
        rfu_rows = rows.loc[rows["rfu_label"].eq(rfu)]
        sequence_ids = set(rfu_unique["unique_sequence_id"])
        evidence_matched_sequences = sequence_ids.intersection(evidence_ids)
        labelled_sequences = sequence_ids.intersection(labelled_ids)
        contributing_sequences = sequence_ids.intersection(contributing_ids)
        matched_cells = rfu_rows.loc[rfu_rows["unique_sequence_id"].isin(evidence_ids)]
        cell_identifier = "cell_id" if "cell_id" in rfu_rows else "input_row_id"
        total_cells = int(rfu_rows[cell_identifier].nunique())
        matched_cell_count = int(matched_cells[cell_identifier].nunique())
        rfu_abundance = abundance.loc[abundance["rfu_label"].eq(rfu)]
        values = (
            rfu_abundance.set_index(antigen_key)["antigen_abundance"]
            if len(rfu_abundance)
            else pd.Series(dtype=float)
        )
        total_labelled = float(values.sum())
        entropy = shannon(values)
        richness = len(values)
        normalized_entropy = (
            entropy / math.log(richness) if richness > 1 else 0.0 if richness == 1 else float("nan")
        )
        eligible = len(contributing_sequences) >= min_matched_sequences and total_labelled > 0
        dominant = (
            sorted(values.index[values.eq(values.max())].astype(str))[0] if len(values) else pd.NA
        )
        score = pd.to_numeric(
            pairs.loc[pairs["unique_sequence_id"].isin(sequence_ids), "evidence_score"],
            errors="coerce",
        )
        ambiguous = status.loc[status["unique_sequence_id"].isin(labelled_sequences), "ambiguous"]
        row: dict[str, Any] = {
            "rfu_label": rfu,
            "assignment_policy": assignment_policy,
            "weighting": weighting,
            "ambiguity_policy": ambiguity_policy,
            "total_rfu_sequences": len(rfu_unique),
            "total_rfu_cells": total_cells,
            "vdjdb_matched_sequences": len(evidence_matched_sequences),
            "antigen_labelled_sequences": len(labelled_sequences),
            "coherence_contributing_sequences": len(contributing_sequences),
            "vdjdb_matched_cells": matched_cell_count,
            "sequence_match_rate": len(evidence_matched_sequences) / len(rfu_unique)
            if len(rfu_unique)
            else 0.0,
            "cell_match_rate": matched_cell_count / total_cells if total_cells else 0.0,
            "evidence_record_count": len(pairs.loc[pairs["unique_sequence_id"].isin(sequence_ids)]),
            "antigen_richness": richness,
            "dominant_antigen": dominant if eligible else pd.NA,
            "dominant_antigen_fraction": float(values.max() / total_labelled)
            if eligible
            else float("nan"),
            "antigen_entropy": entropy if eligible else float("nan"),
            "normalized_antigen_entropy": normalized_entropy if eligible else float("nan"),
            "antigen_purity": float(values.max() / total_labelled) if eligible else float("nan"),
            "ambiguous_sequence_fraction": float(ambiguous.mean()) if len(ambiguous) else 0.0,
            "mean_evidence_score": float(score.mean()) if score.notna().any() else float("nan"),
            "max_evidence_score": float(score.max()) if score.notna().any() else float("nan"),
            "eligible_for_coherence": eligible,
        }
        for output_name, key in (
            ("represented_samples", sample_key),
            ("represented_donors", donor_key),
        ):
            if key is not None:
                if key not in rfu_rows:
                    raise ValueError(f"RFU results are missing metadata column {key!r}.")
                row[output_name] = int(rfu_rows[key].dropna().nunique())
        output.append(row)
    return pd.DataFrame(output)


def _same_antigen_pair_fraction(
    rfu_results: Any,
    evidence: pd.DataFrame,
    *,
    assignment_policy: str,
    antigen_key: str,
    ambiguity_policy: str,
) -> float:
    _, unique = _assigned_unique(rfu_results, assignment_policy=assignment_policy)
    labels, _ = _sequence_antigens(
        evidence, antigen_key=antigen_key, ambiguity_policy=ambiguity_policy
    )
    distributions = {
        sequence_id: dict(zip(group[antigen_key].astype(str), group["label_weight"], strict=True))
        for sequence_id, group in labels.groupby("unique_sequence_id", observed=True)
    }
    numerator = 0.0
    denominator = 0
    for _, group in unique.groupby("rfu_label", observed=True):
        sequence_ids = [value for value in group["unique_sequence_id"] if value in distributions]
        for left_index in range(len(sequence_ids)):
            for right_index in range(left_index + 1, len(sequence_ids)):
                left = distributions[sequence_ids[left_index]]
                right = distributions[sequence_ids[right_index]]
                if ambiguity_policy == "multi_label":
                    similarity = float(bool(set(left).intersection(right)))
                else:
                    similarity = sum(
                        left.get(label, 0.0) * right.get(label, 0.0)
                        for label in set(left).intersection(right)
                    )
                numerator += similarity
                denominator += 1
    return numerator / denominator if denominator else float("nan")


def global_antigen_coherence(
    rfu_results: Any,
    evidence: pd.DataFrame,
    *,
    assignment_policy: str = "nearest",
    antigen_key: str = "epitope",
    ambiguity_policy: str = "fractional",
) -> dict[str, Any]:
    """Return dataset-level descriptive RFU/antigen association measures."""
    _, unique = _assigned_unique(rfu_results, assignment_policy=assignment_policy)
    coherence = rfu_antigen_coherence(
        rfu_results,
        evidence,
        assignment_policy=assignment_policy,
        weighting="unique_sequence",
        antigen_key=antigen_key,
        min_matched_sequences=1,
        ambiguity_policy=ambiguity_policy,
    )
    eligible = coherence.loc[coherence["eligible_for_coherence"]]
    weights = eligible["coherence_contributing_sequences"].astype(float)
    weighted_purity = (
        float(np.average(eligible["antigen_purity"], weights=weights))
        if weights.sum()
        else float("nan")
    )
    weighted_entropy = (
        float(np.average(eligible["antigen_entropy"], weights=weights))
        if weights.sum()
        else float("nan")
    )
    labels, _ = _sequence_antigens(
        evidence, antigen_key=antigen_key, ambiguity_policy=ambiguity_policy
    )
    labelled_ids = set(labels["unique_sequence_id"]).intersection(unique["unique_sequence_id"])
    sizes = unique.groupby("rfu_label", observed=True)["unique_sequence_id"].size()
    multi_ids = set(unique.loc[unique["rfu_label"].map(sizes).ge(2), "unique_sequence_id"])
    fraction_multi = (
        len(labelled_ids.intersection(multi_ids)) / len(labelled_ids)
        if labelled_ids
        else float("nan")
    )
    rfu_antigens = (
        unique[["unique_sequence_id", "rfu_label"]]
        .merge(labels, on="unique_sequence_id", how="inner")
        .groupby("rfu_label", observed=True)[antigen_key]
        .agg(lambda values: set(values.astype(str)))
    )
    shared = 0
    total_pairs = 0
    sets = rfu_antigens.tolist()
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            shared += bool(sets[i].intersection(sets[j]))
            total_pairs += 1
    sharing = shared / total_pairs if total_pairs else float("nan")
    joint = unique[["unique_sequence_id", "rfu_label"]].merge(
        labels, on="unique_sequence_id", how="inner"
    )
    if joint.empty:
        mi = nmi = float("nan")
    else:
        table = joint.pivot_table(
            index="rfu_label",
            columns=antigen_key,
            values="label_weight",
            aggfunc="sum",
            fill_value=0.0,
        )
        probabilities = table.to_numpy(dtype=float)
        probabilities /= probabilities.sum()
        row_p = probabilities.sum(axis=1, keepdims=True)
        col_p = probabilities.sum(axis=0, keepdims=True)
        expected = row_p @ col_p
        positive = probabilities > 0
        mi = float(
            (probabilities[positive] * np.log(probabilities[positive] / expected[positive])).sum()
        )
        h_rfu = float(-(row_p[row_p > 0] * np.log(row_p[row_p > 0])).sum())
        h_antigen = float(-(col_p[col_p > 0] * np.log(col_p[col_p > 0])).sum())
        denominator = math.sqrt(h_rfu * h_antigen)
        nmi = mi / denominator if denominator else float("nan")
    return {
        "weighted_mean_rfu_antigen_purity": weighted_purity,
        "weighted_mean_rfu_antigen_entropy": weighted_entropy,
        "fraction_matched_sequences_in_multisequence_rfus": fraction_multi,
        "same_antigen_pair_fraction_within_rfus": _same_antigen_pair_fraction(
            rfu_results,
            evidence,
            assignment_policy=assignment_policy,
            antigen_key=antigen_key,
            ambiguity_policy=ambiguity_policy,
        ),
        "between_rfu_antigen_sharing": sharing,
        "mutual_information": mi,
        "normalized_mutual_information": nmi,
        "matched_sequence_count": len(labelled_ids),
        "rfu_count": int(unique["rfu_label"].nunique()),
        "assignment_policy": assignment_policy,
        "ambiguity_policy": ambiguity_policy,
        "antigen_key": antigen_key,
    }


def rfu_antigen_permutation_test(
    rfu_results: Any,
    evidence: pd.DataFrame,
    *,
    n_permutations: int = 1000,
    random_state: int = 0,
    stratify_by: str | Sequence[str] | None = None,
    metric: str = "same_antigen_pair_fraction",
    assignment_policy: str = "nearest",
    ambiguity_policy: str = "fractional",
    antigen_key: str = "epitope",
) -> AntigenPermutationResult:
    """Benchmark coherence against size-preserving RFU-label permutations."""
    if (
        isinstance(n_permutations, bool)
        or not isinstance(n_permutations, int)
        or n_permutations < 1
    ):
        raise ValueError("n_permutations must be a positive integer.")
    if metric not in {"same_antigen_pair_fraction", "weighted_mean_purity"}:
        raise ValueError("Unsupported permutation metric.")
    rows, unique = _assigned_unique(rfu_results, assignment_policy=assignment_policy)
    labelled = set(evidence.loc[evidence[antigen_key].notna(), "unique_sequence_id"])
    if len(labelled.intersection(unique["unique_sequence_id"])) < 2:
        raise ValueError("At least two matched unique sequences are required for permutation.")
    fields = (
        []
        if stratify_by is None
        else [stratify_by]
        if isinstance(stratify_by, str)
        else list(stratify_by)
    )
    working = unique.copy()
    normalized_fields: list[str] = []
    for field in fields:
        field_key = str(field).strip().lower()
        if field_key in {"cdr3_length", "length"}:
            working["cdr3_length"] = working["_cdr3_norm"].str.len()
            normalized_fields.append("cdr3_length")
        elif field_key in {"trbv", "v_call"}:
            working["trbv"] = working["_v_norm"].fillna("<missing>")
            normalized_fields.append("trbv")
        elif field not in working:
            raise ValueError(f"Unknown permutation stratification field {field!r}.")
        else:
            normalized_fields.append(str(field))

    def value(frame: pd.DataFrame) -> float:
        if metric == "same_antigen_pair_fraction":
            return _same_antigen_pair_fraction(
                frame,
                evidence,
                assignment_policy="nearest",
                antigen_key=antigen_key,
                ambiguity_policy=ambiguity_policy,
            )
        return global_antigen_coherence(
            frame,
            evidence,
            assignment_policy="nearest",
            antigen_key=antigen_key,
            ambiguity_policy=ambiguity_policy,
        )["weighted_mean_rfu_antigen_purity"]

    observed = value(working)
    if not np.isfinite(observed):
        raise ValueError("Observed permutation metric is undefined for the supplied groups.")
    rng = np.random.default_rng(random_state)
    original_sizes = working["rfu_label"].value_counts().sort_index()
    null = np.empty(n_permutations, dtype=float)
    for permutation in range(n_permutations):
        permuted = working.copy()
        if normalized_fields:
            for _, indices in working.groupby(
                normalized_fields, dropna=False, sort=True
            ).groups.items():
                positions = list(indices)
                permuted.loc[positions, "rfu_label"] = rng.permutation(
                    working.loc[positions, "rfu_label"].to_numpy()
                )
        else:
            permuted["rfu_label"] = rng.permutation(working["rfu_label"].to_numpy())
        if not permuted["rfu_label"].value_counts().sort_index().equals(original_sizes):
            raise RuntimeError("Permutation failed to preserve RFU group sizes.")
        null[permutation] = value(permuted)
    finite = null[np.isfinite(null)]
    if len(finite) != len(null):
        raise ValueError("Permutation metric became undefined under the null model.")
    null_mean = float(finite.mean())
    null_std = float(finite.std(ddof=0))
    probability = (1 + int((finite >= observed).sum())) / (1 + len(finite))
    z_score = (observed - null_mean) / null_std if null_std else float("nan")
    return AntigenPermutationResult(
        observed,
        null,
        null_mean,
        null_std,
        probability,
        z_score,
        {
            "metric": metric,
            "n_permutations": n_permutations,
            "random_state": random_state,
            "stratify_by": normalized_fields,
            "eligible_sequence_count": len(working),
            "matched_sequence_count": len(labelled.intersection(working["unique_sequence_id"])),
            "rfu_count": int(working["rfu_label"].nunique()),
            "assignment_policy": assignment_policy,
            "ambiguity_policy": ambiguity_policy,
            "finite_sample_correction": "(1 + null >= observed) / (1 + n_permutations)",
            "group_sizes_preserved": True,
        },
    )


def compare_antigen_groupings(
    receptor_results: Any,
    evidence: pd.DataFrame,
    *,
    groupings: Sequence[str] = ("rfu", "trbv", "cdr3_length"),
    metrics: Sequence[str] = ("purity", "entropy", "same_antigen_pair_fraction"),
    random_state: int = 0,
    assignment_policy: str = "nearest",
    ambiguity_policy: str = "fractional",
    antigen_key: str = "epitope",
) -> pd.DataFrame:
    """Compare RFU coherence with simple, explicitly limited receptor baselines."""
    _, unique = _assigned_unique(receptor_results, assignment_policy=assignment_policy)
    rng = np.random.default_rng(random_state)
    rows: list[dict[str, Any]] = []
    for grouping in groupings:
        grouped = unique.copy()
        if grouping == "rfu":
            pass
        elif grouping == "trbv":
            grouped["rfu_label"] = grouped["_v_norm"].fillna("<missing>")
        elif grouping == "cdr3_length":
            grouped["rfu_label"] = (
                grouped["_cdr3_norm"].str.len().map(lambda value: f"length_{value}")
            )
        elif grouping in {"trbv_cdr3_length", "trbv+length"}:
            grouped["rfu_label"] = (
                grouped["_v_norm"].fillna("<missing>")
                + "|"
                + grouped["_cdr3_norm"].str.len().astype(str)
            )
        elif grouping in {"random", "size_matched_random"}:
            grouped["rfu_label"] = rng.permutation(grouped["rfu_label"].to_numpy())
        else:
            raise ValueError(f"Unknown antigen grouping baseline {grouping!r}.")
        global_metrics = global_antigen_coherence(
            grouped,
            evidence,
            assignment_policy="nearest",
            antigen_key=antigen_key,
            ambiguity_policy=ambiguity_policy,
        )
        values = {
            "purity": global_metrics["weighted_mean_rfu_antigen_purity"],
            "entropy": global_metrics["weighted_mean_rfu_antigen_entropy"],
            "same_antigen_pair_fraction": global_metrics["same_antigen_pair_fraction_within_rfus"],
        }
        for metric in metrics:
            if metric not in values:
                raise ValueError(f"Unknown grouping comparison metric {metric!r}.")
            rows.append(
                {
                    "grouping_method": grouping,
                    "group_count": int(grouped["rfu_label"].nunique()),
                    "matched_sequence_count": global_metrics["matched_sequence_count"],
                    "metric": metric,
                    "value": values[metric],
                    "weighting": "unique_sequence",
                    "assignment_policy": assignment_policy,
                    "ambiguity_policy": ambiguity_policy,
                    "random_state": random_state,
                }
            )
    return pd.DataFrame(rows)


def summarize_antigen_context(
    rfu_results: Any,
    evidence: pd.DataFrame,
    *,
    metadata: pd.DataFrame | None = None,
    cell_key: str = "cell_id",
    sample_key: str | None = None,
    donor_key: str | None = None,
    phenotype_key: str | None = None,
    condition_key: str | None = None,
    antigen_key: str = "epitope",
) -> AntigenContextResult:
    """Summarize external evidence across explicitly joined biological metadata."""
    rows = analysis_frame(rfu_results).copy()
    if cell_key not in rows:
        raise ValueError(f"RFU results are missing cell identifier {cell_key!r}.")
    before = len(rows)
    unmatched = 0
    if metadata is not None:
        if cell_key not in metadata:
            raise ValueError(f"Metadata are missing cell identifier {cell_key!r}.")
        if metadata[cell_key].isna().any() or metadata[cell_key].duplicated().any():
            raise ValueError("Metadata cell identifiers must be unique and non-missing.")
        overlap = [column for column in metadata if column != cell_key and column in rows]
        if overlap:
            raise ValueError(f"Metadata columns conflict with RFU results: {overlap}")
        unmatched = len(set(rows[cell_key].astype(str)).difference(metadata[cell_key].astype(str)))
        rows = rows.merge(metadata, on=cell_key, how="left", sort=False, validate="many_to_one")
    if len(rows) != before:
        raise RuntimeError("Metadata join expanded antigen-context rows.")
    requested = [key for key in (sample_key, donor_key, phenotype_key, condition_key) if key]
    missing = [key for key in requested if key not in rows]
    if missing:
        raise ValueError(f"Antigen context is missing metadata columns: {missing}")
    if sample_key:
        for key in (donor_key, condition_key):
            if (
                key
                and (rows.groupby(sample_key, observed=True)[key].nunique(dropna=True) > 1).any()
            ):
                raise ValueError(
                    f"Conflicting {key!r} labels occur within at least one {sample_key!r}."
                )
    query = _prepare_queries(rows, chain=None, v_gene_mode="strip_allele")
    matches = evidence.drop_duplicates(["unique_sequence_id", "reference_row_id"])
    evidence_rows = query.merge(
        matches,
        on="unique_sequence_id",
        how="inner",
        suffixes=("", "_evidence"),
        validate="many_to_many",
    )
    sequence_prevalence = pd.DataFrame()
    if sample_key:
        sequence_prevalence = (
            evidence_rows.groupby("unique_sequence_id", observed=True)[sample_key]
            .nunique()
            .rename("sample_prevalence_count")
            .reset_index()
        )
    recurrence_group = [column for column in ("rfu_label", antigen_key) if column in evidence_rows]
    recurrence = pd.DataFrame()
    if recurrence_group:
        aggregations: dict[str, tuple[str, str]] = {}
        if sample_key:
            aggregations["sample_recurrence_count"] = (sample_key, "nunique")
        if donor_key:
            aggregations["donor_recurrence_count"] = (donor_key, "nunique")
        if aggregations:
            recurrence = (
                evidence_rows.groupby(recurrence_group, observed=True)
                .agg(**aggregations)
                .reset_index()
            )
    phenotype_abundance = phenotype_tiers = phenotype_ambiguity = pd.DataFrame()
    if phenotype_key:
        phenotype_abundance = (
            evidence_rows.groupby([phenotype_key, "rfu_label", antigen_key], observed=True)[
                cell_key
            ]
            .nunique()
            .rename("matched_cell_count")
            .reset_index()
        )
        phenotype_tiers = (
            evidence_rows.groupby([phenotype_key, "evidence_tier"], observed=True)[cell_key]
            .nunique()
            .rename("matched_cell_count")
            .reset_index()
        )
        ambiguity = summarize_vdjdb_evidence(rows, evidence).row_summary[
            ["input_row_id", "has_vdjdb_evidence", "antigen_ambiguity"]
        ]
        ambiguity = ambiguity.loc[ambiguity["has_vdjdb_evidence"]].copy()
        phenotype_ambiguity = (
            query[["input_row_id", phenotype_key]]
            .merge(ambiguity, on="input_row_id", validate="one_to_one")
            .groupby(phenotype_key, observed=True)["antigen_ambiguity"]
            .agg(["count", "mean"])
            .reset_index()
        )
    return AntigenContextResult(
        sequence_prevalence,
        recurrence,
        phenotype_abundance,
        phenotype_tiers,
        phenotype_ambiguity,
        {
            "input_row_count": before,
            "output_row_count": len(rows),
            "unmatched_metadata_rows": unmatched,
        },
    )


__all__ = [
    "AntigenContextResult",
    "AntigenPermutationResult",
    "VDJDB_REFERENCE_SCHEMA_VERSION",
    "VDJdbEvidenceSummary",
    "VDJdbReference",
    "annotate_vdjdb",
    "compare_antigen_groupings",
    "global_antigen_coherence",
    "load_vdjdb_reference",
    "normalize_vdjdb_cdr3",
    "normalize_vdjdb_v_gene",
    "rfu_antigen_abundance",
    "rfu_antigen_coherence",
    "rfu_antigen_permutation_test",
    "summarize_antigen_context",
    "summarize_vdjdb_evidence",
    "validate_vdjdb_reference",
]
