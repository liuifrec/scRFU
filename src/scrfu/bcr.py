"""Experimental BCR preprocessing and interpretable state features.

This module does not implement or invoke TCR RFU assignment and does not define
a validated BCR functional-unit reference.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from .pp import normalize_chain, normalize_productive, normalize_text

BCR_SCHEMA_VERSION = "0.1-experimental"
BCR_CHAINS = frozenset({"IGH", "IGK", "IGL"})
BCR_CANONICAL_FIELDS = (
    "cell_id",
    "sequence_id",
    "chain",
    "cdr3aa",
    "cdr3nt",
    "v_call",
    "d_call",
    "j_call",
    "c_call",
    "isotype",
    "productive",
    "paired_receptor_id",
    "clonotype_id",
    "clonal_family_id",
    "mutation_count",
    "mutation_frequency",
    "germline_identity",
    "umi_count",
    "read_count",
    "source_row_id",
)

_ALIASES: dict[str, tuple[str, ...]] = {
    "cell_id": ("cell_id", "cell", "barcode", "cell_barcode"),
    "sequence_id": ("sequence_id", "contig_id", "sequence", "rearrangement_id"),
    "chain": ("chain", "locus", "gene"),
    "cdr3aa": ("cdr3aa", "junction_aa", "cdr3_aa", "cdr3"),
    "cdr3nt": ("cdr3nt", "junction", "junction_nt", "cdr3_nt"),
    "v_call": ("v_call", "v_gene", "v"),
    "d_call": ("d_call", "d_gene", "d"),
    "j_call": ("j_call", "j_gene", "j"),
    "c_call": ("c_call", "c_gene", "constant", "constant_region"),
    "isotype": ("isotype", "subclass", "constant_isotype"),
    "productive": ("productive", "is_productive"),
    "paired_receptor_id": ("paired_receptor_id", "pair_id", "receptor_id"),
    "clonotype_id": ("clonotype_id", "clone_id", "raw_clonotype_id"),
    "clonal_family_id": ("clonal_family_id", "clonal_family", "family_id"),
    "mutation_count": ("mutation_count", "mutations", "v_mutation_count"),
    "mutation_frequency": (
        "mutation_frequency",
        "mutation_freq",
        "shm_frequency",
        "v_mutation_frequency",
    ),
    "germline_identity": ("germline_identity", "v_identity", "v_identity_fraction"),
    "umi_count": ("umi_count", "umis", "umi"),
    "read_count": ("read_count", "reads", "read"),
    "source_row_id": ("source_row_id", "row_id", "source_id"),
}


@dataclass(frozen=True)
class BCRPreparationResult:
    receptors: pd.DataFrame
    pairs: pd.DataFrame
    qc: dict[str, Any]
    provenance: dict[str, Any]


@dataclass(frozen=True)
class BCRFeatureMatrixResult:
    """Experimental one-row-per-cell BCR feature matrix plus missingness QC."""

    features: pd.DataFrame
    missingness: pd.DataFrame
    qc: dict[str, Any]
    parameters: dict[str, Any]


def _header(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def normalize_bcr_chain(value: Any) -> Any:
    """Normalize a BCR chain and reject TCR/unknown loci."""
    chain = normalize_chain(value)
    return chain if pd.isna(chain) or chain in BCR_CHAINS else pd.NA


def normalize_isotype(value: Any) -> Any:
    """Normalize explicit constant-region/isotype text conservatively."""
    normalized = normalize_text(value)
    if pd.isna(normalized):
        return pd.NA
    text = str(normalized).upper().replace(" ", "").replace("-", "")
    text = re.sub(r"\*[A-Z0-9]+$", "", text)
    text = text.removeprefix("IGH")
    text = text.removeprefix("IG")
    match = re.fullmatch(r"(M|D|G[1-4]?|A[12]?|E)", text)
    return f"Ig{match.group(1).title()}" if match else pd.NA


def _fraction(value: Any) -> Any:
    if pd.isna(value):
        return pd.NA
    text = str(value).strip()
    percent = text.endswith("%")
    if percent:
        text = text[:-1]
    try:
        number = float(text)
    except ValueError:
        return pd.NA
    if percent:
        number /= 100.0
    return number if np.isfinite(number) and 0 <= number <= 1 else pd.NA


def _resolve_columns(
    frame: pd.DataFrame, column_mapping: Mapping[str, str] | None
) -> dict[str, str | None]:
    mapping = dict(column_mapping or {})
    unknown = sorted(set(mapping).difference(BCR_CANONICAL_FIELDS))
    if unknown:
        raise ValueError(f"Unknown canonical BCR fields in column_mapping: {unknown}")
    normalized = {_header(column): str(column) for column in frame.columns}
    resolved: dict[str, str | None] = {}
    for field in BCR_CANONICAL_FIELDS:
        if field in mapping:
            if mapping[field] not in frame:
                raise ValueError(f"BCR source column {mapping[field]!r} is missing.")
            resolved[field] = mapping[field]
        else:
            resolved[field] = next(
                (
                    normalized[_header(alias)]
                    for alias in _ALIASES[field]
                    if _header(alias) in normalized
                ),
                None,
            )
    return resolved


def canonicalize_bcr_table(
    data: pd.DataFrame,
    *,
    column_mapping: Mapping[str, str] | None = None,
    source_label: str = "bcr_table",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Return an experimental canonical BCR receptor table and field provenance."""
    if not isinstance(data, pd.DataFrame):
        raise TypeError("BCR input must be a pandas DataFrame.")
    resolved = _resolve_columns(data, column_mapping)
    missing_required = [
        field for field in ("cell_id", "chain", "cdr3aa") if resolved[field] is None
    ]
    if missing_required:
        raise ValueError(f"BCR input is missing required fields: {missing_required}")
    output = pd.DataFrame(index=data.index)
    field_provenance: dict[str, dict[str, Any]] = {}
    for field in BCR_CANONICAL_FIELDS:
        source = resolved[field]
        if source is None:
            output[field] = pd.NA
            field_provenance[field] = {"status": "unavailable", "source_field": None}
        else:
            output[field] = data[source]
            field_provenance[field] = {"status": "direct", "source_field": source}
    if resolved["source_row_id"] is None:
        output["source_row_id"] = data.index.astype(str)
        field_provenance["source_row_id"] = {
            "status": "derived",
            "source_field": "dataframe_index",
        }
    if resolved["sequence_id"] is None:
        output["sequence_id"] = [f"bcr_sequence_{index:09d}" for index in range(len(output))]
        field_provenance["sequence_id"] = {
            "status": "derived",
            "source_field": "stable_source_order",
        }
    for field in (
        "cell_id",
        "sequence_id",
        "cdr3aa",
        "cdr3nt",
        "v_call",
        "d_call",
        "j_call",
        "c_call",
        "paired_receptor_id",
        "clonotype_id",
        "clonal_family_id",
        "source_row_id",
    ):
        output[field] = output[field].map(normalize_text).astype("string")
    output["chain"] = output["chain"].map(normalize_bcr_chain).astype("string")
    output["productive"] = pd.array(output["productive"].map(normalize_productive), dtype="boolean")
    explicit_isotype = output["isotype"].map(normalize_isotype)
    constant_isotype = output["c_call"].map(normalize_isotype)
    output["isotype"] = explicit_isotype.fillna(constant_isotype).astype("string")
    if resolved["isotype"] is None and output["isotype"].notna().any():
        field_provenance["isotype"] = {"status": "derived", "source_field": resolved["c_call"]}
    output["mutation_count"] = pd.to_numeric(output["mutation_count"], errors="coerce").astype(
        "Float64"
    )
    output["umi_count"] = pd.to_numeric(output["umi_count"], errors="coerce").astype("Float64")
    output["read_count"] = pd.to_numeric(output["read_count"], errors="coerce").astype("Float64")
    output["mutation_frequency"] = output["mutation_frequency"].map(_fraction).astype("Float64")
    output["germline_identity"] = output["germline_identity"].map(_fraction).astype("Float64")
    derive_mutation = output["mutation_frequency"].isna() & output["germline_identity"].notna()
    output.loc[derive_mutation, "mutation_frequency"] = (
        1.0 - output.loc[derive_mutation, "germline_identity"]
    )
    if derive_mutation.any() and resolved["mutation_frequency"] is None:
        field_provenance["mutation_frequency"] = {
            "status": "derived",
            "source_field": resolved["germline_identity"],
            "rule": "1 - germline_identity",
        }
    output.insert(0, "input_row_id", [f"bcr_row_{index:09d}" for index in range(len(output))])
    output.insert(1, "source_order", np.arange(len(output), dtype=int))
    if output["cell_id"].isna().any() or output["sequence_id"].isna().any():
        raise ValueError("Canonical BCR cell_id and sequence_id values must be non-missing.")
    if output["sequence_id"].duplicated().any():
        raise ValueError("Canonical BCR sequence_id values must be unique.")
    provenance = {
        "schema_version": BCR_SCHEMA_VERSION,
        "experimental": True,
        "source_label": source_label,
        "field_provenance": field_provenance,
        "tcr_rfu_assignment_permitted": False,
    }
    return output.reset_index(drop=True), provenance


def select_productive_bcr_chains(table: pd.DataFrame) -> pd.DataFrame:
    """Mark deterministic primary heavy and light records while retaining every row."""
    frame = table.copy().reset_index(drop=True)
    required = {"cell_id", "sequence_id", "chain", "productive", "source_order"}
    if not required.issubset(frame):
        raise ValueError(f"BCR chain selection requires fields: {sorted(required)}")
    frame["bcr_arm"] = (
        frame["chain"].map({"IGH": "heavy", "IGK": "light", "IGL": "light"}).astype("string")
    )
    frame["selected_for_pair"] = False
    frame["selection_rank"] = pd.Series(pd.NA, index=frame.index, dtype="Int64")
    productive_rank = frame["productive"].map({True: 0, False: 2}).fillna(1).astype(int)
    umi = pd.to_numeric(
        frame.get("umi_count", pd.Series(np.nan, index=frame.index)), errors="coerce"
    )
    reads = pd.to_numeric(
        frame.get("read_count", pd.Series(np.nan, index=frame.index)), errors="coerce"
    )
    ranked = frame.assign(
        _productive_rank=productive_rank,
        _umi_rank=umi.fillna(-1),
        _read_rank=reads.fillna(-1),
    ).sort_values(
        ["cell_id", "bcr_arm", "_productive_rank", "_umi_rank", "_read_rank", "source_order"],
        ascending=[True, True, True, False, False, True],
        kind="stable",
    )
    valid = ranked["bcr_arm"].notna()
    ranked.loc[valid, "selection_rank"] = (
        ranked.loc[valid].groupby(["cell_id", "bcr_arm"], observed=True).cumcount() + 1
    )
    selected_indices = ranked.loc[valid & ranked["selection_rank"].eq(1)].index
    frame.loc[selected_indices, "selected_for_pair"] = True
    frame.loc[ranked.index, "selection_rank"] = ranked["selection_rank"]
    return frame


def pair_bcr_chains(table: pd.DataFrame) -> pd.DataFrame:
    """Return one deterministic heavy/light pairing record per cell."""
    selected = select_productive_bcr_chains(table)
    chosen = selected.loc[selected["selected_for_pair"]].copy()
    heavy = chosen.loc[chosen["bcr_arm"].eq("heavy")].set_index("cell_id")
    light = chosen.loc[chosen["bcr_arm"].eq("light")].set_index("cell_id")
    cells = pd.Index(sorted(selected["cell_id"].dropna().astype(str).unique()), name="cell_id")
    pairs = pd.DataFrame(index=cells)
    pairs["heavy_sequence_id"] = heavy["sequence_id"].reindex(cells)
    pairs["heavy_chain"] = heavy["chain"].reindex(cells)
    pairs["light_sequence_id"] = light["sequence_id"].reindex(cells)
    pairs["light_chain"] = light["chain"].reindex(cells)
    pairs["pair_status"] = np.select(
        [
            pairs["heavy_sequence_id"].notna() & pairs["light_sequence_id"].notna(),
            pairs["heavy_sequence_id"].notna(),
            pairs["light_sequence_id"].notna(),
        ],
        ["paired", "heavy_only", "light_only"],
        default="no_bcr_chain",
    )
    pairs["paired_receptor_id"] = pd.Series(pd.NA, index=pairs.index, dtype="string")
    paired = pairs["pair_status"].eq("paired")
    pairs.loc[paired, "paired_receptor_id"] = "bcr_pair:" + pairs.index[paired].astype(str)
    counts = selected.groupby(["cell_id", "bcr_arm"], observed=True).size().unstack(fill_value=0)
    pairs["heavy_candidate_count"] = counts.get("heavy", pd.Series(0, index=cells)).reindex(
        cells, fill_value=0
    )
    pairs["light_candidate_count"] = counts.get("light", pd.Series(0, index=cells)).reindex(
        cells, fill_value=0
    )
    return pairs.reset_index()


def bcr_state_features(table: pd.DataFrame, pairs: pd.DataFrame | None = None) -> pd.DataFrame:
    """Add conservative receptor-derived BCR state features without phenotype inference."""
    frame = table.copy().reset_index(drop=True)
    required = {"cell_id", "chain", "cdr3aa", "isotype", "clonotype_id", "clonal_family_id"}
    if not required.issubset(frame):
        raise ValueError(f"BCR state features require fields: {sorted(required)}")
    frame["cdr3_length"] = frame["cdr3aa"].astype("string").str.len().astype("Int64")
    frame["class_switched"] = pd.Series(pd.NA, index=frame.index, dtype="boolean")
    frame.loc[frame["isotype"].isin(["IgM", "IgD"]), "class_switched"] = False
    frame.loc[
        frame["isotype"].astype("string").str.match(r"^Ig(?:G[1-4]?|A[12]?|E)$", na=False),
        "class_switched",
    ] = True
    pairing = pair_bcr_chains(frame) if pairs is None else pairs.copy()
    frame = frame.merge(
        pairing[["cell_id", "pair_status", "light_chain"]],
        on="cell_id",
        how="left",
        sort=False,
        validate="many_to_one",
    )
    frame["paired_light_chain_status"] = frame["pair_status"].eq("paired")
    for key, output in (
        ("clonotype_id", "clonotype_size"),
        ("clonal_family_id", "clonal_family_size"),
    ):
        sizes = frame.loc[frame[key].notna()].groupby(key, observed=True)["sequence_id"].nunique()
        frame[output] = frame[key].map(sizes).astype("Int64")
    family = frame.loc[frame["clonal_family_id"].notna()].groupby("clonal_family_id", observed=True)
    richness = family["cdr3aa"].nunique()
    size = family["sequence_id"].nunique()
    frame["within_family_cdr3_richness"] = frame["clonal_family_id"].map(richness).astype("Int64")
    frame["within_family_cdr3_diversity"] = (
        frame["clonal_family_id"].map(richness) / frame["clonal_family_id"].map(size)
    ).astype("Float64")
    return frame


def bcr_feature_matrix(
    table: pd.DataFrame,
    *,
    pairs: pd.DataFrame | None = None,
) -> BCRFeatureMatrixResult:
    """Build an interpretable experimental BCR feature matrix per cell.

    Features are receptor-derived only. Missing SHM, isotype, light-chain, or
    clonal-family values remain explicitly missing; no phenotype or outcome
    label is inferred or used to fit this representation.
    """
    pairing = pair_bcr_chains(table) if pairs is None else pairs.copy()
    states = bcr_state_features(table, pairing)
    selected = states.loc[states["selected_for_pair"]].copy()
    common = [
        "cell_id",
        "source_row_id",
        "cdr3aa",
        "cdr3_length",
        "v_call",
        "j_call",
        "isotype",
        "class_switched",
        "mutation_count",
        "mutation_frequency",
        "germline_identity",
        "clonotype_id",
        "clonotype_size",
        "clonal_family_id",
        "clonal_family_size",
        "within_family_cdr3_diversity",
    ]

    def arm_frame(arm: str) -> pd.DataFrame:
        arm_rows = selected.loc[selected["bcr_arm"].eq(arm), common].copy()
        if arm_rows["cell_id"].duplicated().any():
            raise RuntimeError(
                f"Deterministic BCR selection produced multiple {arm} rows per cell."
            )
        return arm_rows.rename(
            columns={column: f"{arm}_{column}" for column in common if column != "cell_id"}
        )

    features = pairing.copy()
    features = features.merge(arm_frame("heavy"), on="cell_id", how="left", validate="one_to_one")
    features = features.merge(arm_frame("light"), on="cell_id", how="left", validate="one_to_one")
    if len(features) != pairing["cell_id"].nunique() or features["cell_id"].duplicated().any():
        raise RuntimeError("BCR feature construction changed the one-row-per-cell contract.")
    feature_columns = [
        column
        for column in features
        if column
        not in {
            "cell_id",
            "pair_status",
            "paired_receptor_id",
            "heavy_candidate_count",
            "light_candidate_count",
        }
    ]
    missingness = pd.DataFrame(
        {
            "feature": feature_columns,
            "available_count": [int(features[column].notna().sum()) for column in feature_columns],
            "missing_count": [int(features[column].isna().sum()) for column in feature_columns],
        }
    )
    missingness["available_fraction"] = missingness["available_count"] / max(1, len(features))
    qc = {
        "schema_version": BCR_SCHEMA_VERSION,
        "experimental": True,
        "cell_count": len(features),
        "feature_count": len(feature_columns),
        "paired_cell_count": int(features["pair_status"].eq("paired").sum()),
        "heavy_only_cell_count": int(features["pair_status"].eq("heavy_only").sum()),
        "light_only_cell_count": int(features["pair_status"].eq("light_only").sum()),
        "tcr_rfu_assignment_permitted": False,
    }
    parameters = {
        "unit": "cell",
        "heavy_chain": "IGH",
        "light_chains": ["IGK", "IGL"],
        "selection": "productive, UMI, read count, source order",
        "outcome_labels_used": False,
        "functional_unit_reference": None,
    }
    return BCRFeatureMatrixResult(features.reset_index(drop=True), missingness, qc, parameters)


def bcr_qc_summary(table: pd.DataFrame, pairs: pd.DataFrame | None = None) -> dict[str, Any]:
    """Summarize experimental BCR chain, pairing, and field completeness QC."""
    pairing = pair_bcr_chains(table) if pairs is None else pairs
    return {
        "schema_version": BCR_SCHEMA_VERSION,
        "experimental": True,
        "receptor_rows": len(table),
        "cell_count": int(table["cell_id"].nunique()),
        "chain_counts": {
            str(key): int(value) for key, value in table["chain"].value_counts(dropna=False).items()
        },
        "productive_rows": int(table["productive"].eq(True).sum()),
        "isotype_available_rows": int(table["isotype"].notna().sum()),
        "mutation_frequency_available_rows": int(table["mutation_frequency"].notna().sum()),
        "clonal_family_available_rows": int(table["clonal_family_id"].notna().sum()),
        "pair_status_counts": {
            str(key): int(value) for key, value in pairing["pair_status"].value_counts().items()
        },
        "tcr_rfu_assignment_permitted": False,
    }


def prepare_bcr_table(
    data: pd.DataFrame,
    *,
    column_mapping: Mapping[str, str] | None = None,
    source_label: str = "bcr_table",
) -> BCRPreparationResult:
    """Canonicalize, retain, pair, and summarize BCR receptor records.

    This experimental function never invokes TCR RFU assignment and does not
    construct a BCR functional reference.
    """
    receptors, provenance = canonicalize_bcr_table(
        data, column_mapping=column_mapping, source_label=source_label
    )
    selected = select_productive_bcr_chains(receptors)
    pairs = pair_bcr_chains(selected)
    return BCRPreparationResult(selected, pairs, bcr_qc_summary(selected, pairs), provenance)


__all__ = [
    "BCR_CANONICAL_FIELDS",
    "BCR_CHAINS",
    "BCR_SCHEMA_VERSION",
    "BCRPreparationResult",
    "BCRFeatureMatrixResult",
    "bcr_feature_matrix",
    "bcr_qc_summary",
    "bcr_state_features",
    "canonicalize_bcr_table",
    "normalize_bcr_chain",
    "normalize_isotype",
    "pair_bcr_chains",
    "prepare_bcr_table",
    "select_productive_bcr_chains",
]
