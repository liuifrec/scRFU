from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .pp import canonicalize_receptor_table, normalize_productive, normalize_text

ADAPTER_API_VERSION = "1"


@dataclass(frozen=True)
class AdapterResult:
    receptors: pd.DataFrame
    cell_metadata: pd.DataFrame
    qc: dict[str, Any]
    adapter_name: str
    adapter_version: str
    provenance: dict[str, Any]


ReceptorAdapter = Callable[..., AdapterResult]

_CELLRANGER_COLUMNS: dict[str, str] = {
    "barcode": "cell_id",
    "chain": "chain",
    "cdr3": "cdr3aa",
    "v_gene": "v_call",
    "d_gene": "d_call",
    "j_gene": "j_call",
    "c_gene": "c_call",
    "cdr3_nt": "junction",
    "umis": "umi_count",
    "reads": "read_count",
    "raw_clonotype_id": "clonotype_id",
}

_AIRR_ALIASES: dict[str, tuple[str, ...]] = {
    "cell_id": ("cell_id", "cell", "barcode", "cellid"),
    "chain": ("chain", "locus"),
    "cdr3aa": ("cdr3aa", "junction_aa", "cdr3_aa", "cdr3"),
    "v_call": ("v_call", "v_gene", "trbv", "v"),
    "productive": ("productive",),
    "duplicate_count": ("duplicate_count", "duplicates"),
    "consensus_count": ("consensus_count", "consensus"),
    "umi_count": ("umi_count", "umis", "umi"),
    "read_count": ("read_count", "reads", "junction_reads"),
    "j_call": ("j_call", "j_gene"),
    "d_call": ("d_call", "d_gene"),
    "c_call": ("c_call", "c_gene"),
    "junction": ("junction",),
    "junction_aa": ("junction_aa",),
    "clonotype_id": ("clonotype_id", "clone_id"),
}


def _pick_columns(frame: pd.DataFrame, aliases: Mapping[str, Sequence[str]]) -> dict[str, Any]:
    normalized = {str(column).strip().lower(): column for column in frame.columns}
    return {
        target: next((normalized[name] for name in names if name in normalized), None)
        for target, names in aliases.items()
    }


def _metadata_frame(
    source: Any,
    *,
    metadata: pd.DataFrame | None,
    metadata_columns: Sequence[str],
) -> pd.DataFrame:
    if metadata is not None:
        result = metadata.copy()
    elif hasattr(source, "obs") and isinstance(source.obs, pd.DataFrame):
        result = source.obs.copy()
    else:
        return pd.DataFrame(index=pd.Index([], name="cell_id"))
    result.index = pd.Index(result.index.astype(str), name="cell_id")
    if result.index.has_duplicates or result.index.isna().any():
        raise ValueError("Cell metadata must have unique, non-missing cell identifiers.")
    missing = [column for column in metadata_columns if column not in result]
    if missing:
        raise ValueError(f"Requested metadata columns are missing: {missing}")
    return result.loc[:, list(metadata_columns)] if metadata_columns else result.iloc[:, :0]


def _source_frame(source: Any, *, airr_key: str) -> tuple[pd.DataFrame, str]:
    if isinstance(source, pd.DataFrame):
        return source.copy(), "dataframe"
    if hasattr(source, "obsm") and airr_key in source.obsm:
        value = source.obsm[airr_key]
        if isinstance(value, pd.DataFrame):
            return value.copy(), f"obsm[{airr_key!r}]"
        if value.__class__.__module__.startswith("awkward"):
            raise TypeError(
                "Awkward-array scirpy AIRR objects are not yet supported directly; "
                "convert the receptor records to a pandas DataFrame first."
            )
        try:
            return pd.DataFrame(value), f"obsm[{airr_key!r}]"
        except Exception as exc:
            raise ValueError(f"Could not convert obsm[{airr_key!r}] to a DataFrame.") from exc
    raise KeyError(f"Input has no DataFrame source and no obsm[{airr_key!r}] receptor table.")


def _numeric_rank(values: pd.Series) -> pd.Series:
    return pd.to_numeric(values, errors="coerce").fillna(float("-inf"))


def adapt_airr_dataframe(
    source: Any,
    *,
    adapter_name: str = "generic_airr_dataframe",
    airr_key: str = "airr",
    chain: str | None = None,
    productive_only: bool = True,
    primary_chain: bool = True,
    metadata: pd.DataFrame | None = None,
    metadata_columns: Sequence[str] = (),
) -> AdapterResult:
    """Adapt a pandas/scirpy AIRR-like table to canonical receptor rows."""
    frame, source_format = _source_frame(source, airr_key=airr_key)
    columns = _pick_columns(frame, _AIRR_ALIASES)
    missing = [name for name in ("cell_id", "chain", "cdr3aa", "v_call") if columns[name] is None]
    if missing:
        raise ValueError(
            f"AIRR table missing columns: {missing}. Columns seen: {list(frame.columns)}"
        )

    selected = pd.DataFrame(index=frame.index)
    selected["cell_id"] = frame[columns["cell_id"]]
    selected["chain"] = frame[columns["chain"]]
    selected["cdr3aa"] = frame[columns["cdr3aa"]]
    selected["v_call"] = frame[columns["v_call"]]
    selected["productive"] = (
        frame[columns["productive"]]
        if columns["productive"] is not None
        else pd.Series(pd.NA, index=frame.index, dtype="boolean")
    )
    for optional in _AIRR_ALIASES:
        if optional in selected or columns[optional] is None:
            continue
        selected[optional] = frame[columns[optional]]
    selected["source_row_id"] = pd.Series(frame.index.astype(str), index=frame.index)
    selected["_source_order"] = range(len(selected))
    selected = canonicalize_receptor_table(selected, source_adapter=adapter_name)

    requested_chain = str(chain).strip().upper() if chain is not None else None
    if requested_chain is not None:
        selected = selected.loc[selected["chain"].eq(requested_chain)].copy()
    if productive_only and columns["productive"] is not None:
        selected = selected.loc[selected["productive"].fillna(False)].copy()

    cell_metadata = _metadata_frame(source, metadata=metadata, metadata_columns=metadata_columns)
    if len(cell_metadata.index) == 0:
        cell_metadata = pd.DataFrame(
            index=pd.Index(
                selected["cell_id"].dropna().astype(str).drop_duplicates(), name="cell_id"
            )
        )
    if not cell_metadata.empty or len(cell_metadata.index):
        receptor_cells = set(selected["cell_id"].dropna().astype(str))
        unknown = sorted(receptor_cells.difference(cell_metadata.index.astype(str)))
        if unknown:
            raise ValueError(
                "AIRR receptor barcodes are absent from cell metadata/obs: "
                + ", ".join(unknown[:10])
            )

    before_primary = len(selected)
    if primary_chain and not selected.empty:
        rank_columns: list[str] = []
        selected["_rank_productive"] = selected["productive"].fillna(False).astype(int)
        rank_columns.append("_rank_productive")
        for column in ("umi_count", "read_count", "duplicate_count", "consensus_count"):
            if column in selected:
                rank = f"_rank_{column}"
                selected[rank] = _numeric_rank(selected[column])
                rank_columns.append(rank)
        selected = selected.sort_values(
            ["cell_id", *rank_columns, "_source_order"],
            ascending=[True, *([False] * len(rank_columns)), True],
            kind="stable",
        ).drop_duplicates("cell_id", keep="first")
        selected = selected.sort_values("_source_order", kind="stable")
        selected = selected.drop(columns=rank_columns)

    selected = selected.drop(columns="_source_order").reset_index(drop=True)
    selected["input_row_id"] = [f"row_{index:08d}" for index in range(len(selected))]
    selected = canonicalize_receptor_table(selected)
    qc = {
        "source_row_count": len(frame),
        "selected_receptor_row_count": len(selected),
        "selected_cell_count": int(selected["cell_id"].nunique()),
        "rows_before_primary_selection": before_primary,
        "chain_counts": selected["chain"].value_counts().to_dict(),
        "productive_row_count": int(selected["productive"].fillna(False).sum()),
        "primary_chain": primary_chain,
        "productive_only": productive_only,
        "requested_chain": requested_chain,
    }
    return AdapterResult(
        receptors=selected,
        cell_metadata=cell_metadata,
        qc=qc,
        adapter_name=adapter_name,
        adapter_version=ADAPTER_API_VERSION,
        provenance={"source_format": source_format, "airr_key": airr_key},
    )


def _wells_column(frame: pd.DataFrame, slot: int, field: str) -> Any | None:
    normalized = {str(column).lower(): column for column in frame.columns}
    aliases = [f"tcr-ir_vdj_{slot}_{field}", f"ir_vdj_{slot}_{field}"]
    if field == "junction_aa":
        aliases.extend([f"tcr-ir_vdj_{slot}_cdr3", f"ir_vdj_{slot}_cdr3"])
    return next((normalized[alias] for alias in aliases if alias in normalized), None)


def _wells_source(
    source: Any, key: str, metadata_columns: Sequence[str]
) -> tuple[pd.DataFrame, pd.DataFrame, str, dict[str, Any]]:
    if isinstance(source, (str, Path)):
        from .wells import read_wells_receptors_h5ad

        data = read_wells_receptors_h5ad(source, obs_columns=metadata_columns)
        return (
            data.tcr_ir.copy(),
            data.obs.copy(),
            data.tcr_ir_container,
            {
                "source_path": str(Path(source).expanduser().resolve()),
                "source_atlas_dimensions": list(data.atlas_shape),
            },
        )
    if hasattr(source, "tcr_ir") and hasattr(source, "obs"):
        return (
            source.tcr_ir.copy(),
            source.obs.copy(),
            source.tcr_ir_container,
            {"source_atlas_dimensions": list(source.atlas_shape)},
        )
    if hasattr(source, "uns") and key in source.uns:
        value, container = source.uns[key], "uns"
    elif hasattr(source, "obsm") and key in source.obsm:
        value, container = source.obsm[key], "obsm"
    else:
        raise KeyError(f"Input contains neither uns[{key!r}] nor obsm[{key!r}].")
    frame = value.copy() if isinstance(value, pd.DataFrame) else pd.DataFrame(value)
    if not hasattr(source, "obs_names"):
        raise ValueError("Wells input requires obs_names for exact row alignment.")
    obs_names = pd.Index(map(str, source.obs_names), name="cell_id")
    if len(frame) != len(obs_names):
        raise ValueError("Wells TCR_IR row count does not match observations.")
    if isinstance(frame.index, pd.RangeIndex):
        frame.index = obs_names
    elif list(map(str, frame.index)) != obs_names.tolist():
        raise ValueError("Wells TCR_IR index does not exactly match obs_names.")
    metadata = source.obs.copy() if hasattr(source, "obs") else pd.DataFrame(index=obs_names)
    metadata.index = obs_names
    return frame, metadata, container, {}


def adapt_wells_tcr_ir(
    source: Any,
    *,
    key: str = "TCR_IR",
    chain: str | None = "TRB",
    productive_only: bool = True,
    primary_chain: bool = True,
    metadata_columns: Sequence[str] = (),
) -> AdapterResult:
    """Interpret Wells flattened VDJ_1/VDJ_2 fields behind a named adapter."""
    frame, metadata, container, source_provenance = _wells_source(source, key, metadata_columns)
    metadata.index = pd.Index(metadata.index.astype(str), name="cell_id")
    missing_metadata = [column for column in metadata_columns if column not in metadata]
    if missing_metadata:
        raise ValueError(f"Requested metadata columns are missing: {missing_metadata}")
    metadata = metadata.loc[:, list(metadata_columns)] if metadata_columns else metadata.iloc[:, :0]

    rows: list[dict[str, Any]] = []
    for slot in (1, 2):
        cdr3_col = _wells_column(frame, slot, "junction_aa")
        v_col = _wells_column(frame, slot, "v_call")
        if cdr3_col is None or v_col is None:
            if slot == 1:
                raise ValueError("Wells TCR_IR is missing VDJ_1 CDR3 or V-call fields.")
            continue
        for source_order, (cell_id, values) in enumerate(frame.iterrows()):
            locus_col = _wells_column(frame, slot, "locus")
            productive_col = _wells_column(frame, slot, "productive")
            duplicate_col = _wells_column(frame, slot, "duplicate_count")
            consensus_col = _wells_column(frame, slot, "consensus_count")
            cdr3 = normalize_text(values[cdr3_col])
            v_call = normalize_text(values[v_col])
            if pd.isna(cdr3) or pd.isna(v_call):
                continue
            locus = values[locus_col] if locus_col is not None else chain
            productive = (
                normalize_productive(values[productive_col])
                if productive_col is not None
                else pd.NA
            )
            rows.append(
                {
                    "cell_id": str(cell_id),
                    "chain": locus,
                    "cdr3aa": cdr3,
                    "v_call": v_call,
                    "productive": productive,
                    "source_adapter": "wells_tcr_ir",
                    "source_row_id": str(cell_id),
                    "source_slot": f"VDJ_{slot}",
                    "duplicate_count": values[duplicate_col]
                    if duplicate_col is not None
                    else pd.NA,
                    "consensus_count": values[consensus_col]
                    if consensus_col is not None
                    else pd.NA,
                    "_source_order": source_order,
                    "_slot_order": slot,
                }
            )
    selected = canonicalize_receptor_table(pd.DataFrame(rows), source_adapter="wells_tcr_ir")
    requested_chain = str(chain).upper() if chain is not None else None
    if requested_chain is not None and not selected.empty:
        selected = selected.loc[selected["chain"].eq(requested_chain)].copy()
    if productive_only and not selected.empty and selected["productive"].notna().any():
        selected = selected.loc[selected["productive"].fillna(False)].copy()
    before_primary = len(selected)
    if primary_chain and not selected.empty:
        selected["_rank_duplicate"] = _numeric_rank(selected["duplicate_count"])
        selected["_rank_consensus"] = _numeric_rank(selected["consensus_count"])
        selected = selected.sort_values(
            ["cell_id", "_rank_duplicate", "_rank_consensus", "_slot_order"],
            ascending=[True, False, False, True],
            kind="stable",
        ).drop_duplicates("cell_id", keep="first")
        selected = selected.sort_values(["_source_order", "_slot_order"], kind="stable")
        selected = selected.drop(columns=["_rank_duplicate", "_rank_consensus"])
    selected = selected.drop(columns=["_source_order", "_slot_order"]).reset_index(drop=True)
    selected["input_row_id"] = [f"row_{index:08d}" for index in range(len(selected))]
    selected = canonicalize_receptor_table(selected)
    qc = {
        "tcr_ir_source": container,
        "source_row_count": len(frame),
        "selected_receptor_row_count": len(selected),
        "selected_cell_count": int(selected["cell_id"].nunique()) if not selected.empty else 0,
        "rows_before_primary_selection": before_primary,
        "chain_counts": selected["chain"].value_counts().to_dict() if not selected.empty else {},
        "productive_row_count": int(selected["productive"].fillna(False).sum())
        if not selected.empty
        else 0,
        "non_c_rows": int((~selected["cdr3aa"].str.startswith("C", na=False)).sum())
        if not selected.empty
        else 0,
        "primary_chain": primary_chain,
        "productive_only": productive_only,
        "requested_chain": requested_chain,
    }
    return AdapterResult(
        receptors=selected,
        cell_metadata=metadata,
        qc=qc,
        adapter_name="wells_tcr_ir",
        adapter_version=ADAPTER_API_VERSION,
        provenance={"source_format": "h5ad_wells_tcr_ir", "tcr_key": key, **source_provenance},
    )


def _truth_mask(values: pd.Series, *, column: str) -> pd.Series:
    normalized = values.map(normalize_productive).astype("boolean")
    invalid = values.notna() & normalized.isna()
    if invalid.any():
        examples = values.loc[invalid].astype(str).drop_duplicates().tolist()[:5]
        raise ValueError(f"Cell Ranger column {column!r} contains invalid booleans: {examples}")
    return normalized.fillna(False)


def adapt_cellranger_vdj(
    source: pd.DataFrame | str | Path,
    *,
    chain: str | None = None,
    productive_only: bool = True,
    primary_chain: bool = True,
    filter_is_cell: bool = True,
    filter_high_confidence: bool = True,
    metadata: pd.DataFrame | None = None,
    metadata_columns: Sequence[str] = (),
) -> AdapterResult:
    """Adapt a Cell Ranger V(D)J contig annotation table.

    ``source`` may be an in-memory DataFrame or a comma/tab-delimited file.
    Optional Cell Ranger columns are retained when present, while barcode,
    chain, and amino-acid CDR3 are the only mandatory source fields.
    """
    if isinstance(source, pd.DataFrame):
        frame = source.copy()
        source_format = "dataframe"
        source_path: str | None = None
    elif isinstance(source, (str, Path)):
        path = Path(source).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Cell Ranger VDJ table not found: {path}")
        separator = "\t" if ".tsv" in path.name.lower() else ","
        frame = pd.read_csv(path, sep=separator)
        source_format = "cellranger_vdj_csv" if separator == "," else "cellranger_vdj_tsv"
        source_path = str(path)
    else:
        raise TypeError("Cell Ranger input must be a pandas DataFrame or CSV/TSV path.")

    normalized = {str(column).strip().lower(): column for column in frame.columns}
    missing = [name for name in ("barcode", "chain", "cdr3") if name not in normalized]
    if missing:
        raise ValueError(
            "Cell Ranger VDJ table missing required columns "
            f"{missing}; columns seen: {list(frame.columns)}"
        )

    selected = pd.DataFrame(index=frame.index)
    for source_name, target_name in _CELLRANGER_COLUMNS.items():
        if source_name in normalized:
            selected[target_name] = frame[normalized[source_name]]
    if "v_call" not in selected:
        selected["v_call"] = pd.Series(pd.NA, index=frame.index, dtype="string")
    selected["productive"] = (
        frame[normalized["productive"]]
        if "productive" in normalized
        else pd.Series(pd.NA, index=frame.index, dtype="boolean")
    )
    for column in (
        "contig_id",
        "is_cell",
        "high_confidence",
        "length",
        "full_length",
        "raw_consensus_id",
    ):
        if column in normalized:
            selected[column] = frame[normalized[column]]
    selected["source_row_id"] = pd.Series(frame.index.astype(str), index=frame.index)
    selected["_source_order"] = range(len(selected))
    selected = canonicalize_receptor_table(selected, source_adapter="cellranger_vdj")

    source_rows = len(selected)
    selected = selected.loc[
        selected["cell_id"].notna() & selected["chain"].notna() & selected["cdr3aa"].notna()
    ].copy()
    after_required_values = len(selected)
    if filter_is_cell and "is_cell" in selected:
        selected = selected.loc[_truth_mask(selected["is_cell"], column="is_cell")].copy()
    after_is_cell = len(selected)
    if filter_high_confidence and "high_confidence" in selected:
        selected = selected.loc[
            _truth_mask(selected["high_confidence"], column="high_confidence")
        ].copy()
    after_high_confidence = len(selected)
    requested_chain = str(chain).strip().upper() if chain is not None else None
    if requested_chain is not None:
        selected = selected.loc[selected["chain"].eq(requested_chain)].copy()
    if productive_only and "productive" in normalized:
        selected = selected.loc[selected["productive"].fillna(False)].copy()
    before_primary = len(selected)

    if primary_chain and not selected.empty:
        selected["_rank_productive"] = selected["productive"].fillna(False).astype(int)
        rank_columns = ["_rank_productive"]
        if "high_confidence" in selected:
            selected["_rank_confidence"] = _truth_mask(
                selected["high_confidence"], column="high_confidence"
            ).astype(int)
            rank_columns.append("_rank_confidence")
        for column in ("umi_count", "read_count"):
            if column in selected:
                rank = f"_rank_{column}"
                selected[rank] = _numeric_rank(selected[column])
                rank_columns.append(rank)
        selected = selected.sort_values(
            ["cell_id", *rank_columns, "_source_order"],
            ascending=[True, *([False] * len(rank_columns)), True],
            kind="stable",
        ).drop_duplicates("cell_id", keep="first")
        selected = selected.sort_values("_source_order", kind="stable").drop(columns=rank_columns)

    selected = selected.drop(columns="_source_order").reset_index(drop=True)
    selected["input_row_id"] = [f"row_{index:08d}" for index in range(len(selected))]
    selected = canonicalize_receptor_table(selected)
    if metadata is not None:
        cell_metadata = _metadata_frame(frame, metadata=metadata, metadata_columns=metadata_columns)
    elif metadata_columns:
        missing_metadata = [column for column in metadata_columns if column not in frame]
        if missing_metadata:
            raise ValueError(f"Requested metadata columns are missing: {missing_metadata}")
        metadata_rows: list[dict[str, Any]] = []
        for barcode in selected["cell_id"].dropna().astype(str).drop_duplicates():
            source_rows_for_cell = frame.loc[frame[normalized["barcode"]].astype(str).eq(barcode)]
            row: dict[str, Any] = {"cell_id": barcode}
            for column in metadata_columns:
                values = source_rows_for_cell[column].dropna().drop_duplicates().tolist()
                if len(values) > 1:
                    raise ValueError(
                        f"Cell Ranger barcode {barcode!r} has conflicting metadata "
                        f"for {column!r}: {values[:5]}"
                    )
                row[column] = values[0] if values else pd.NA
            metadata_rows.append(row)
        cell_metadata = pd.DataFrame(
            metadata_rows, columns=["cell_id", *metadata_columns]
        ).set_index("cell_id")
    else:
        cell_metadata = pd.DataFrame(
            index=pd.Index(
                selected["cell_id"].dropna().astype(str).drop_duplicates(), name="cell_id"
            )
        )
    unknown = sorted(
        set(selected["cell_id"].dropna().astype(str)).difference(cell_metadata.index.astype(str))
    )
    if unknown:
        raise ValueError(
            "Cell Ranger receptor barcodes are absent from cell metadata: "
            + ", ".join(unknown[:10])
        )

    qc = {
        "source_row_count": source_rows,
        "rows_after_required_value_filter": after_required_values,
        "rows_after_is_cell_filter": after_is_cell,
        "rows_after_high_confidence_filter": after_high_confidence,
        "rows_before_primary_selection": before_primary,
        "selected_receptor_row_count": len(selected),
        "selected_cell_count": int(selected["cell_id"].nunique()),
        "chain_counts": selected["chain"].value_counts().to_dict(),
        "productive_row_count": int(selected["productive"].fillna(False).sum()),
        "filter_is_cell": filter_is_cell,
        "filter_high_confidence": filter_high_confidence,
        "productive_only": productive_only,
        "primary_chain": primary_chain,
        "requested_chain": requested_chain,
    }
    provenance: dict[str, Any] = {"source_format": source_format}
    if source_path is not None:
        provenance["source_path"] = source_path
    return AdapterResult(
        receptors=selected,
        cell_metadata=cell_metadata,
        qc=qc,
        adapter_name="cellranger_vdj",
        adapter_version=ADAPTER_API_VERSION,
        provenance=provenance,
    )


def _scirpy_airr(source: Any, **kwargs: Any) -> AdapterResult:
    return adapt_airr_dataframe(source, adapter_name="scirpy_airr", **kwargs)


def _generic_airr(source: Any, **kwargs: Any) -> AdapterResult:
    return adapt_airr_dataframe(source, adapter_name="generic_airr_dataframe", **kwargs)


_ADAPTERS: dict[str, ReceptorAdapter] = {
    "cellranger_vdj": adapt_cellranger_vdj,
    "wells_tcr_ir": adapt_wells_tcr_ir,
    "scirpy_airr": _scirpy_airr,
    "generic_airr_dataframe": _generic_airr,
}
_ALIASES = {
    "cellranger": "cellranger_vdj",
    "tenx_vdj": "cellranger_vdj",
    "wells": "wells_tcr_ir",
    "airr": "generic_airr_dataframe",
    "dataframe": "generic_airr_dataframe",
}


def get_receptor_adapter(name: str) -> ReceptorAdapter:
    normalized = str(name).strip().lower()
    normalized = _ALIASES.get(normalized, normalized)
    if normalized not in _ADAPTERS:
        raise ValueError(
            f"Unknown receptor adapter {name!r}. Available: {', '.join(list_receptor_adapters())}"
        )
    return _ADAPTERS[normalized]


def list_receptor_adapters() -> list[str]:
    return sorted(_ADAPTERS)


def prepare_receptors(source: Any, *, adapter: str, **kwargs: Any) -> AdapterResult:
    return get_receptor_adapter(adapter)(source, **kwargs)


__all__ = [
    "ADAPTER_API_VERSION",
    "AdapterResult",
    "ReceptorAdapter",
    "adapt_airr_dataframe",
    "adapt_cellranger_vdj",
    "adapt_wells_tcr_ir",
    "get_receptor_adapter",
    "list_receptor_adapters",
    "prepare_receptors",
]
