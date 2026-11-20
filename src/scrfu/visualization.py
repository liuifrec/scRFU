"""Inspectible preparation for descriptive receptor/RFU plots; no plotting imports.

These functions never assign RFUs, call clonotypes, fit an embedding or mutate
their inputs. Result tables retain the observation unit and denominators.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import pandas as pd

from .adapters import AdapterResult, prepare_receptors
from .downstream import _ordered, _require_columns, _sort_key
from .pp import normalize_chain, normalize_productive
from .repertoire import analysis_frame
from .scverse import _airr_records, _resolve_airr_modality, validate_scrfu_schema

PlotUnit = Literal["cell", "receptor", "read"]
_WEIGHT = "__plot_weight"
_IDENTITY = "__plot_identity"
_GROUP = "__plot_group"


@dataclass(frozen=True)
class ReceptorPlotData:
    """Validated rows, count unit, identity definition and preparation audit.

    ``rows`` retains literal calls and source row order after explicit filters.
    Cell weights are one per selected cell; read weights use the supplied count
    column. Unique-receptor weights are deduplicated *within each plotted group*.
    DataFrames remain caller-inspectible; preparation makes defensive copies.
    """

    rows: pd.DataFrame
    unit: PlotUnit
    identity_cols: tuple[str, ...]
    cell_col: str
    audit: dict[str, Any]


@dataclass(frozen=True)
class EmbeddingPlotData:
    """Existing coordinates and exactly aligned metadata, in coordinate order."""

    coordinates: np.ndarray
    metadata: pd.DataFrame
    audit: dict[str, Any]


@dataclass(frozen=True)
class VJUsageResult:
    """Complete literal pairs, displayed marginals and complete group denominators."""

    literal_pairs: pd.DataFrame
    pairs: pd.DataFrame
    v_counts: pd.DataFrame
    j_counts: pd.DataFrame
    denominators: pd.Series
    unit: str
    audit: dict[str, Any]


@dataclass(frozen=True)
class SpectratypeResult:
    """Length-bin counts, valid-sequence denominators and explicit exclusions."""

    counts: pd.DataFrame
    bin_edges: np.ndarray
    denominators: pd.Series
    input_denominators: pd.Series
    excluded: pd.DataFrame
    unit: str
    audit: dict[str, Any]


@dataclass(frozen=True)
class MembershipResult:
    """Complete identity/RFU nodes and confirmed edges; layout is descriptive."""

    receptors: pd.DataFrame
    rfus: pd.DataFrame
    edges: pd.DataFrame
    unit: str
    identity_cols: tuple[str, ...]
    audit: dict[str, Any]


@dataclass(frozen=True)
class LongitudinalPresenceResult:
    """Per-donor identity matrices with unavailable samples represented by NaN.

    Denominators count included observations at each sample, in the selected
    unit. ``presence_flows`` is defined only for two observed samples and counts
    distinct identities in their detected union, independently of that unit.
    """

    matrices: dict[Any, pd.DataFrame]
    denominators: pd.DataFrame
    observed: pd.DataFrame
    presence_flows: pd.DataFrame
    sample_order: tuple[Any, ...]
    identity_cols: tuple[str, ...]
    unit: str
    min_count: float
    audit: dict[str, Any]


def _index(values: Any, *, name: str | None = None) -> pd.Index:
    return pd.Index(list(values), name=name, tupleize_cols=False)


def _ids(values: Any, *, context: str, unique: bool = True) -> pd.Index:
    index = _index(values, name=getattr(values, "name", None))
    if index.isna().any() or any(isinstance(x, str) and not x.strip() for x in index):
        raise ValueError(f"{context} requires non-missing, non-blank identifiers.")
    if unique and index.has_duplicates:
        raise ValueError(f"{context} requires unique identifiers; duplicate cells/records found.")
    return index


def _frame(frame: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(frame, pd.DataFrame):
        raise TypeError("Expected a pandas DataFrame.")
    if not frame.columns.is_unique:
        raise ValueError("Input has duplicate column names.")
    if frame.empty:
        raise ValueError("Plot input is empty.")
    if any(str(c).startswith("__plot_") for c in frame.columns):
        raise ValueError("Column names starting '__plot_' are reserved for plotting preparation.")
    return frame.copy()


def _bool(values: pd.Series, name: str) -> pd.Series:
    result = values.map(normalize_productive)
    invalid = values.notna() & result.isna()
    if invalid.any():
        raise ValueError(f"{name!r} contains invalid boolean values.")
    return result.astype("boolean")


def _join_metadata(rows: pd.DataFrame, metadata: pd.DataFrame, cell_col: str) -> pd.DataFrame:
    if not metadata.columns.is_unique:
        raise ValueError("Metadata has duplicate column names.")
    meta = metadata.copy()
    if cell_col in meta:
        meta = meta.set_index(cell_col)
    _ids(meta.index, context="Cell metadata")
    _require_columns(rows, [cell_col], "Metadata join")
    if not rows[cell_col].isin(meta.index).all():
        raise ValueError(
            "Receptor cell identifiers are absent from metadata; no join was inferred."
        )
    # Never replace chain-level calls/assignments with a cell-level summary.
    for col in meta.columns.difference(rows.columns, sort=False):
        rows[col] = rows[cell_col].map(meta[col])
    return rows


def _source_rows(
    data: Any,
    *,
    airr_mod: str,
    airr_key: str | None,
    assignment_key: str | None,
    cell_col: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if isinstance(data, AdapterResult):
        return _join_metadata(data.receptors.copy(), data.cell_metadata, cell_col), {
            "source": "AdapterResult",
            "adapter": data.adapter_name,
        }
    if isinstance(data, pd.DataFrame):
        return data.copy(), {"source": "table"}
    if airr_key is None:
        if hasattr(data, "mod"):
            data, _ = _resolve_airr_modality(data, airr_mod)
        _ids(data.obs.index, context="Observations")
        rows = analysis_frame(data)
        if cell_col not in rows:
            rows[cell_col] = data.obs.index
        return rows, {"source": "obs", "airr_key": None}
    target, _ = _resolve_airr_modality(data, airr_mod)
    _ids(target.obs.index, context="AIRR observations")
    if airr_key not in target.obsm:
        # Explicit obs-only path; never guess a study-specific receptor field.
        raise KeyError(f"No obsm[{airr_key!r}]; pass airr_key=None to use explicit .obs columns.")
    if isinstance(target.obsm[airr_key], pd.DataFrame):
        rows = target.obsm[airr_key].copy()
        return _join_metadata(rows, target.obs, cell_col), {"source": "AIRR DataFrame"}
    _, chains = _airr_records(target, airr_key)
    assigned = None
    reference = None
    if assignment_key is not None and assignment_key in target.obsm:
        check = validate_scrfu_schema(
            data, airr_mod=airr_mod, airr_key=airr_key, key_added=assignment_key
        )
        reference = check["reference_identity"]
        import awkward as ak

        assigned = ak.to_list(target.obsm[assignment_key])
    # Reuse the canonical adapter without selecting or discarding any chain.
    adapted = prepare_receptors(
        target,
        adapter="scirpy_airr",
        airr_key=airr_key,
        productive_only=False,
        primary_chain=False,
    )
    rows = adapted.receptors.copy()
    raw_records = [record for cell in chains for record in cell]
    if len(raw_records) != len(rows):
        raise ValueError("Canonical AIRR adapter changed chain alignment.")
    for col in sorted({k for record in raw_records for k in record} - set(rows.columns)):
        rows[col] = [record.get(col) for record in raw_records]
    # Reuse canonical identifiers/locus handling while preserving literal calls.
    for dest, fields in {
        "cdr3aa": ("junction_aa", "cdr3aa"),
        "v_call": ("v_call",),
        "j_call": ("j_call",),
        "junction": ("junction",),
        "productive": ("productive",),
    }.items():
        if any(any(key in record for key in fields) for record in raw_records):
            rows[dest] = [
                next((record[key] for key in fields if key in record), None)
                for record in raw_records
            ]
    if cell_col != "cell_id":
        rows = rows.rename(columns={"cell_id": cell_col})
    if assigned is not None:
        flat = []
        for records, results in zip(chains, assigned, strict=True):
            for i, (record, result) in enumerate(zip(records, results, strict=True)):
                locus = normalize_chain(record.get("locus", record.get("chain")))
                saved_locus = normalize_chain(result.get("locus"))
                same_locus = (pd.isna(locus) and pd.isna(saved_locus)) or (
                    not pd.isna(locus) and not pd.isna(saved_locus) and locus == saved_locus
                )
                if result["chain_index"] != i or not same_locus:
                    raise ValueError("Saved RFU chain indices/loci disagree with the AIRR records.")
                flat.append(result)
        if len(flat) != len(rows):
            raise ValueError("Canonical AIRR adapter changed chain alignment.")
        for col in ("rfu_label", "rfu_score", "threshold_pass", "eligible", "reference_identity"):
            rows[col] = [r[col] for r in flat]
    elif "rfu_label" in target.obs and "rfu_label" not in rows:
        raise ValueError(
            "Only cell-level RFUs are available; do not broadcast them to AIRR chains. "
            "Use airr_key=None for an explicitly selected cell table."
        )
    return _join_metadata(rows, target.obs, cell_col), {
        "source": "native AIRR",
        "reference_identity": reference,
        "chain_alignment_verified": assigned is not None,
    }


def prepare_receptor_plot_data(
    data: Any,
    *,
    unit: PlotUnit = "cell",
    identity_cols: Sequence[str] = (),
    cell_col: str = "cell_id",
    chain: str | None = "TRB",
    chain_col: str = "chain",
    productive_only: bool = False,
    productive_col: str = "productive",
    multi_chain: str = "error",
    primary_col: str | None = None,
    read_count_col: str | None = None,
    record_id_col: str = "input_row_id",
    qualified_only: bool = False,
    qualification_col: str = "pass_thr",
    metadata: pd.DataFrame | None = None,
    airr_mod: str = "airr",
    airr_key: str | None = "airr",
    assignment_key: str | None = "scrfu",
) -> ReceptorPlotData:
    """Prepare canonical, AIRR, AnnData or MuData rows without assigning RFUs.

    Defaults select literal TRB rows but do not assume productivity or score
    qualification. Set ``chain=None`` for a preselected table lacking a chain
    column. Native AIRR calls use the existing canonical adapter and chain-aligned
    scRFU storage; ``airr_key=None`` explicitly uses observation columns.

    ``multi_chain='error'`` rejects duplicate cells after locus/productivity
    filtering. ``'primary'`` requires exactly one saved True primary flag per
    represented cell. ``'all'`` is allowed only for receptor/read weighting.
    Cell identifiers must be globally unique across samples after selection;
    construct an explicit composite identifier when barcodes repeat.

    Read weighting requires finite nonnegative integer counts AND unique record
    IDs; neither UMIs nor row multiplicity are inferred as sequencing reads.
    Receptor weighting requires explicit non-missing identity columns and counts
    each identity once per group. Identity tuples never imply paired clonotypes.
    Metadata joins use identifiers, never row position. Existing receptor columns
    take precedence over cell metadata; inputs are not modified.
    """
    if unit not in {"cell", "receptor", "read"}:
        raise ValueError("unit must be 'cell', 'receptor', or 'read'.")
    if multi_chain not in {"error", "primary", "all"}:
        raise ValueError("multi_chain must be 'error', 'primary', or 'all'.")
    if unit == "cell" and multi_chain == "all":
        raise ValueError("Cell weighting cannot count all chains; supply an explicit primary flag.")
    identity_cols = tuple(identity_cols)
    if unit == "receptor" and not identity_cols:
        raise ValueError("Receptor weighting requires explicit identity_cols.")
    rows, audit = _source_rows(
        data, airr_mod=airr_mod, airr_key=airr_key, assignment_key=assignment_key, cell_col=cell_col
    )
    rows = _frame(rows)
    audit["input_rows"] = len(rows)
    if metadata is not None:
        rows = _join_metadata(rows, metadata, cell_col)
    if chain is not None:
        _require_columns(rows, [chain_col], "Chain filtering (or pass chain=None)")
        rows = rows.loc[
            rows[chain_col].map(normalize_chain).eq(normalize_chain(chain)).fillna(False)
        ].copy()
    audit["after_chain_filter"] = len(rows)
    if productive_only:
        _require_columns(rows, [productive_col], "Productivity filtering")
        rows = rows.loc[_bool(rows[productive_col], productive_col).fillna(False)].copy()
    audit["after_productivity_filter"] = len(rows)
    if rows.empty:
        raise ValueError("No receptor rows remain after chain/productivity filtering.")
    if cell_col in rows:
        _ids(rows[cell_col], context="Receptor cells", unique=False)
        if multi_chain == "primary":
            if primary_col is None:
                raise ValueError("multi_chain='primary' requires primary_col.")
            _require_columns(rows, [primary_col], "Primary selection")
            selected = _bool(rows[primary_col], primary_col).fillna(False)
            n = selected.groupby(rows[cell_col], observed=True, sort=False).sum()
            if not n.eq(1).all():
                raise ValueError("Exactly one saved primary receptor is required per cell.")
            rows = rows.loc[selected].copy()
        elif multi_chain == "error":
            _ids(rows[cell_col], context="Single-chain cell weighting")
    elif unit == "cell" or multi_chain == "primary":
        raise ValueError(f"Cell/primary selection requires {cell_col!r}.")
    if record_id_col in rows:
        _ids(rows[record_id_col], context="Receptor record IDs")
    if unit == "read":
        if read_count_col is None:
            raise ValueError("Read weighting requires read_count_col and unique record IDs.")
        _require_columns(rows, [read_count_col, record_id_col], "Read weighting")
        weights = pd.to_numeric(rows[read_count_col], errors="raise").to_numpy(dtype=float)
        if (
            not np.isfinite(weights).all()
            or (weights < 0).any()
            or not np.equal(weights, np.floor(weights)).all()
        ):
            raise ValueError("Read counts must be finite, nonnegative integers.")
        rows[_WEIGHT] = weights
    else:
        rows[_WEIGHT] = 1.0
    if qualified_only:
        _require_columns(rows, [qualification_col], "Qualification filtering")
        rows = rows.loc[_bool(rows[qualification_col], qualification_col).fillna(False)].copy()
    if rows.empty or not rows[_WEIGHT].gt(0).any():
        raise ValueError("No positive observations remain after the explicit filters.")
    if identity_cols:
        _identity(rows, identity_cols)
    audit.update(
        selected_rows=len(rows),
        selected_cells=int(rows[cell_col].nunique()) if cell_col in rows else None,
        selected_row_weight=float(rows[_WEIGHT].sum()),
        unit=unit,
        identity_cols=list(identity_cols),
        chain=chain,
        productive_only=productive_only,
        multi_chain=multi_chain,
        primary_col=primary_col,
        qualified_only=qualified_only,
        qualification_col=qualification_col if qualified_only else None,
        read_count_col=read_count_col if unit == "read" else None,
        metadata_precedence="receptor columns retained; only additional metadata columns joined",
    )
    return ReceptorPlotData(rows.reset_index(drop=True), unit, identity_cols, cell_col, audit)


def _prepared(
    data: Any, unit: PlotUnit | None, identity_cols: Sequence[str] | None
) -> ReceptorPlotData:
    if isinstance(data, ReceptorPlotData):
        if unit is not None and unit != data.unit:
            raise ValueError(
                "unit disagrees with prepared data; prepare the requested unit explicitly."
            )
        if identity_cols is not None and tuple(identity_cols) != data.identity_cols:
            raise ValueError("identity_cols disagree with prepared data.")
        return data
    return prepare_receptor_plot_data(data, unit=unit or "cell", identity_cols=identity_cols or ())


def _identity(frame: pd.DataFrame, cols: Sequence[str]) -> pd.Series:
    if not cols or len(set(cols)) != len(cols):
        raise ValueError("identity_cols must contain distinct explicit columns.")
    _require_columns(frame, cols, "Observed identity")
    if frame[list(cols)].isna().any().any() or frame[list(cols)].eq("").any().any():
        raise ValueError(
            "Observed identity columns contain missing/blank values; resolve or explicitly filter them."
        )
    return pd.Series(
        list(frame[list(cols)].itertuples(index=False, name=None)), index=frame.index, dtype=object
    )


def _weighted(
    prepared: ReceptorPlotData, frame: pd.DataFrame, groups: list[str], values: list[str]
) -> pd.DataFrame:
    frame = frame.copy()
    if prepared.unit == "receptor":
        frame[_IDENTITY] = _identity(frame, prepared.identity_cols)
        keys = list(dict.fromkeys([*groups, _IDENTITY]))
        diversity = frame.groupby(keys, sort=False, observed=True, dropna=False)[values].nunique(
            dropna=False
        )
        if diversity.gt(1).any().any():
            raise ValueError(
                "One receptor identity has ambiguous plotted values within a group; "
                "use a finer identity definition or cell/read weighting."
            )
        frame = frame.drop_duplicates(keys).copy()
        frame[_WEIGHT] = 1.0
    return frame


def _grouped(
    frame: pd.DataFrame, groupby: str | None, group_order: Sequence[Any] | None
) -> tuple[pd.DataFrame, list[Any]]:
    frame = frame.copy()
    if groupby is None:
        frame[_GROUP] = "All observations"
    else:
        _require_columns(frame, [groupby], "Grouping")
        if frame[groupby].isna().any():
            raise ValueError(
                "Grouping values are missing; label or explicitly filter them before preparation."
            )
        frame[_GROUP] = frame[groupby]
    groups = _ordered(frame[_GROUP])
    if group_order is not None:
        order = _ids(group_order, context="group_order").tolist()
        if set(order) != set(groups):
            raise ValueError("group_order must contain every observed group exactly once.")
        groups = order
    return frame, groups


def prepare_rfu_embedding(
    data: Any = None,
    *,
    coordinates: Any = None,
    metadata: pd.DataFrame | None = None,
    cell_ids: Sequence[Any] | None = None,
    cell_col: str | None = None,
    basis: str = "X_umap",
    modality: str | None = None,
    allow_missing_metadata: bool = False,
) -> EmbeddingPlotData:
    """Align saved 2-D coordinates and cell metadata by exact identifiers.

    AnnData uses its observation order. MuData requires an explicit ``modality``.
    A coordinate DataFrame supplies its own cell index; an ndarray requires
    ``cell_ids``. Extra metadata cells are rejected. Missing metadata cells raise
    unless explicitly allowed, in which case they remain missing and are audited.
    No embedding, sorting, subsampling or expression access occurs.
    """
    if hasattr(data, "mod"):
        if modality is None or modality not in data.mod:
            raise ValueError("Select an existing MuData modality explicitly.")
        data = data.mod[modality]
    if data is not None and hasattr(data, "obs"):
        _ids(data.obs.index, context="Embedding observations")
        if coordinates is None:
            if basis not in data.obsm:
                raise KeyError(
                    f"Saved embedding obsm[{basis!r}] is unavailable; no embedding is computed."
                )
            coordinates = data.obsm[basis]
        if cell_ids is None:
            cell_ids = data.obs.index
        if metadata is None:
            metadata = data.obs
    elif isinstance(data, pd.DataFrame) and metadata is None:
        metadata = data
    if isinstance(coordinates, pd.DataFrame):
        coordinate_ids = _ids(coordinates.index, context="Coordinate rows")
        if cell_ids is not None and not coordinate_ids.equals(_index(cell_ids)):
            raise ValueError("Coordinate DataFrame index and cell_ids disagree in order.")
        cell_ids = coordinate_ids
    if coordinates is None or cell_ids is None or metadata is None:
        raise ValueError(
            "Existing coordinates, explicit cell identifiers and metadata are required."
        )
    ids = _ids(cell_ids, context="Embedding cells")
    xy = np.asarray(coordinates, dtype=float)
    if xy.shape != (len(ids), 2) or len(ids) == 0 or not np.isfinite(xy).all():
        raise ValueError("Coordinates must be a nonempty finite (n_cells, 2) array.")
    if not metadata.columns.is_unique:
        raise ValueError("Metadata has duplicate column names.")
    meta = metadata.copy()
    if cell_col is not None:
        _require_columns(meta, [cell_col], "Embedding metadata")
        meta = meta.set_index(cell_col)
    _ids(meta.index, context="Embedding metadata")
    if not meta.index.isin(ids).all():
        raise ValueError("Metadata contains cell identifiers absent from coordinates.")
    missing = ~ids.isin(meta.index)
    if missing.any() and not allow_missing_metadata:
        raise ValueError("Coordinate cells are missing from metadata; no join was inferred.")
    return EmbeddingPlotData(
        xy.copy(),
        meta.reindex(ids).copy(),
        {
            "source_cells": len(ids),
            "missing_metadata_cells": int(missing.sum()),
            "basis": basis,
            "modality": modality,
            "embedding_computed": False,
            "coordinate_order_preserved": True,
        },
    )


def _calls(frame: pd.DataFrame, col: str, missing_label: str) -> pd.Series:
    values = frame[col].copy() if col in frame else pd.Series(pd.NA, index=frame.index)
    if values.dropna().eq(missing_label).any():
        raise ValueError(f"Literal call collides with reserved category {missing_label!r}.")
    return values.astype(object).mask(values.isna() | values.eq(""), missing_label)


def _categories(
    values: pd.Series, weights: pd.Series, top_n: int | None, missing: str, other: str
) -> tuple[list[Any], dict[Any, Any]]:
    if values.eq(other).any():
        raise ValueError(f"Literal call collides with reserved category {other!r}.")
    calls = _ordered(values[values.ne(missing)])
    if top_n is not None:
        if isinstance(top_n, bool) or not isinstance(top_n, int) or top_n < 1:
            raise ValueError("top_n must be a positive integer or None.")
        totals = weights.groupby(values, observed=True).sum()
        calls = sorted(calls, key=lambda x: (-totals[x], _sort_key(x), str(x)))[:top_n]
        calls = _ordered(calls)
    mapping = {v: v if v in calls or v == missing else other for v in values.unique()}
    categories = calls + ([missing] if values.eq(missing).any() else [])
    if other in mapping.values():
        categories.append(other)
    return categories, mapping


def prepare_rfu_vj_usage(
    data: Any,
    *,
    groupby: str | None = None,
    unit: PlotUnit | None = None,
    identity_cols: Sequence[str] | None = None,
    v_col: str = "v_call",
    j_col: str = "j_call",
    top_n: int | None = None,
    group_order: Sequence[Any] | None = None,
) -> VJUsageResult:
    """Count literal V/J pairs and marginals without discarding Other or missing mass.

    Raw input uses the default single-TRB preparation; supply ReceptorPlotData for
    other locus/primary/qualification/read/metadata policies. No top-N limit is
    applied by default. A limit selects V and J separately by total weight across
    groups, breaking ties naturally then lexically, and aggregates remaining calls
    into explicit Other categories. Missing fields/NA/blank calls remain separate.
    Receptor weighting rejects identities mapping to multiple V/J pairs per group.
    """
    prepared = _prepared(data, unit, identity_cols)
    frame, groups = _grouped(prepared.rows, groupby, group_order)
    if v_col not in frame and j_col not in frame:
        raise ValueError("V/J usage requires at least one source V or J call column.")
    frame["__plot_v"] = _calls(frame, v_col, "Missing V")
    frame["__plot_j"] = _calls(frame, j_col, "Missing J")
    weighted = _weighted(prepared, frame, [_GROUP], ["__plot_v", "__plot_j"])
    den = weighted.groupby(_GROUP, observed=True)[_WEIGHT].sum().reindex(groups)
    v_categories, vm = _categories(
        weighted["__plot_v"], weighted[_WEIGHT], top_n, "Missing V", "Other V"
    )
    j_categories, jm = _categories(
        weighted["__plot_j"], weighted[_WEIGHT], top_n, "Missing J", "Other J"
    )
    literal = (
        weighted.groupby([_GROUP, "__plot_v", "__plot_j"], observed=True, sort=False)[_WEIGHT]
        .sum()
        .reset_index()
    )
    literal.columns = ["group", "v_call", "j_call", "value"]
    display = literal.assign(v_call=literal.v_call.map(vm), j_call=literal.j_call.map(jm))
    display = (
        display.groupby(["group", "v_call", "j_call"], observed=True, sort=False)
        .value.sum()
        .reset_index()
    )

    def matrix(col: str, categories: list[Any]) -> pd.DataFrame:
        return (
            display.groupby(["group", col], observed=True)
            .value.sum()
            .unstack(fill_value=0)
            .reindex(index=groups, columns=categories, fill_value=0)
        )

    v, j = matrix("v_call", v_categories), matrix("j_call", j_categories)
    audit = dict(
        prepared.audit,
        groupby=groupby,
        top_n=top_n,
        denominator="all selected observations in each group, including missing/Other calls",
        missing_columns=[c for c in [v_col, j_col] if c not in prepared.rows],
        v_mapping=vm,
        j_mapping=jm,
    )
    return VJUsageResult(literal, display, v, j, den, prepared.unit, audit)


def prepare_rfu_cdr3_spectratype(
    data: Any,
    *,
    groupby: str | None = None,
    unit: PlotUnit | None = None,
    identity_cols: Sequence[str] | None = None,
    sequence_col: str = "cdr3aa",
    bins: int | Sequence[float] | None = None,
    invalid: str = "error",
    group_order: Sequence[Any] | None = None,
) -> SpectratypeResult:
    """Bin lengths of actual uppercase, unambiguous 20-letter amino-acid sequences.

    Missing/blank sequences are excluded and counted. Invalid sequences raise by
    default; ``invalid='exclude'`` records their excluded weight. No trimming,
    translation or padding is performed. Duplicate identities receive one vote
    per group only in explicit receptor mode. Normalization uses valid sequences;
    input denominators and excluded mass are returned separately. Explicit bin
    edges must cover all valid lengths; the final right edge is inclusive.
    """
    if invalid not in {"error", "exclude"}:
        raise ValueError("invalid must be 'error' or 'exclude'.")
    prepared = _prepared(data, unit, identity_cols)
    frame, groups = _grouped(prepared.rows, groupby, group_order)
    _require_columns(frame, [sequence_col], "Spectratype")
    frame = _weighted(prepared, frame, [_GROUP], [sequence_col])
    missing = frame[sequence_col].isna() | frame[sequence_col].eq("")
    valid = frame[sequence_col].map(
        lambda s: isinstance(s, str) and bool(s) and all(c in "ACDEFGHIKLMNPQRSTVWY" for c in s)
    )
    bad = ~missing & ~valid
    if bad.any() and invalid == "error":
        raise ValueError(
            "Invalid amino-acid sequences; expected uppercase unambiguous amino acids."
        )
    lengths = frame.loc[valid, sequence_col].str.len().to_numpy(dtype=int)
    if len(lengths) == 0:
        raise ValueError("No valid amino-acid sequences remain for the spectratype.")
    if bins is None:
        edges = np.arange(lengths.min() - 0.5, lengths.max() + 1.5)
    else:
        if isinstance(bins, bool) or (isinstance(bins, int) and bins < 1):
            raise ValueError("bins must be positive or strictly increasing explicit edges.")
        edges = np.histogram_bin_edges(lengths, bins=bins)
    if (
        not np.isfinite(edges).all()
        or not np.all(np.diff(edges) > 0)
        or edges[0] > lengths.min()
        or edges[-1] < lengths.max()
    ):
        raise ValueError(
            "Bin edges must be finite, strictly increasing and cover every valid length."
        )
    table = pd.DataFrame(
        0.0, index=_index(groups, name="group"), columns=pd.RangeIndex(len(edges) - 1, name="bin")
    )
    excluded = pd.DataFrame(0.0, index=table.index, columns=["missing", "invalid"])
    for g in groups:
        in_group = frame[_GROUP].eq(g)
        selected = frame.loc[in_group & valid]
        table.loc[g] = np.histogram(
            selected[sequence_col].str.len(), bins=edges, weights=selected[_WEIGHT]
        )[0]
        excluded.loc[g] = [
            frame.loc[in_group & missing, _WEIGHT].sum(),
            frame.loc[in_group & bad, _WEIGHT].sum(),
        ]
    input_den = frame.groupby(_GROUP, observed=True)[_WEIGHT].sum().reindex(groups)
    return SpectratypeResult(
        table,
        edges,
        table.sum(axis=1),
        input_den,
        excluded,
        prepared.unit,
        dict(
            prepared.audit,
            sequence_col=sequence_col,
            invalid=invalid,
            groupby=groupby,
            denominator="valid sequence weight within group; exclusions and input denominator returned separately",
        ),
    )


def prepare_rfu_membership(
    data: Any,
    *,
    identity_cols: Sequence[str] = ("cdr3aa", "v_call"),
    rfu_col: str = "rfu_label",
    unit: PlotUnit | None = None,
    donor_col: str | None = None,
    allow_multiple_rfus: bool = False,
) -> MembershipResult:
    """Prepare confirmed identity-to-RFU edges and complete abundance counts.

    The default identity is literal amino-acid CDR3 plus V call. Supply nucleotide
    and gene columns, or an independently established clonotype column, explicitly
    for another definition. Missing RFUs are excluded and audited; missing identity
    components raise. Conflicting assignments for one identity raise unless
    explicitly allowed. Mixed non-missing reference identities always raise.
    No distances, clustering, lineage or functional equivalence are inferred.
    """
    prepared = _prepared(
        data,
        unit,
        identity_cols if unit == "receptor" and not isinstance(data, ReceptorPlotData) else None,
    )
    if prepared.unit == "receptor" and tuple(identity_cols) != prepared.identity_cols:
        raise ValueError("Membership identity_cols must match receptor-weighted preparation.")
    frame = prepared.rows.copy()
    _require_columns(frame, [rfu_col], "Membership")
    if "reference_identity" in frame and frame.reference_identity.dropna().nunique() > 1:
        raise ValueError("Membership cannot combine incompatible reference identities.")
    missing = frame[rfu_col].isna() | frame[rfu_col].eq("")
    excluded = (
        float(_identity(frame.loc[missing], prepared.identity_cols).nunique())
        if prepared.unit == "receptor"
        else float(frame.loc[missing, _WEIGHT].sum())
    )
    frame = frame.loc[~missing].copy()
    if frame.empty:
        raise ValueError("No assigned RFU memberships remain.")
    frame[_IDENTITY] = _identity(frame, identity_cols)
    if (
        not allow_multiple_rfus
        and frame.groupby(_IDENTITY, sort=False)[rfu_col].nunique().gt(1).any()
    ):
        raise ValueError(
            "An observed identity has multiple RFU assignments; resolve reference/assignment ambiguity."
        )
    if donor_col is not None:
        _require_columns(frame, [donor_col], "Membership donor representation")
        if frame[donor_col].isna().any():
            raise ValueError("Membership donor representation contains missing donors.")
    original = frame.copy()
    weighted = _weighted(prepared, frame, [rfu_col], [rfu_col])
    weighted[_IDENTITY] = _identity(weighted, identity_cols)
    identities = sorted(
        weighted[_IDENTITY].unique(), key=lambda t: tuple((type(x).__name__, str(x)) for x in t)
    )
    rfu_labels = _ordered(weighted[rfu_col])
    ri = {v: f"receptor_{i}" for i, v in enumerate(identities)}
    gi = {v: f"rfu_{i}" for i, v in enumerate(rfu_labels)}
    edges = (
        weighted.groupby([_IDENTITY, rfu_col], observed=True, sort=False)[_WEIGHT]
        .sum()
        .reset_index()
    )
    edges.columns = ["identity", "rfu", "value"]
    edges["receptor_node"] = edges.identity.map(ri)
    edges["rfu_node"] = edges.rfu.map(gi)
    abundance = weighted.groupby(_IDENTITY, sort=False)[_WEIGHT].sum()
    if prepared.unit == "receptor":
        abundance[:] = 1.0
    receptors = pd.DataFrame(
        {
            "node_id": list(ri.values()),
            "identity": identities,
            "abundance": [abundance[x] for x in identities],
        }
    )
    if donor_col is not None:
        donors = original.groupby(_IDENTITY, sort=False)[donor_col].nunique()
        receptors["donors"] = [donors[x] for x in identities]
    rfus = pd.DataFrame({"node_id": list(gi.values()), "rfu": rfu_labels})
    rfus["abundance"] = rfus.rfu.map(edges.groupby("rfu", observed=True).value.sum())
    rfus["identities"] = rfus.rfu.map(edges.groupby("rfu", observed=True).receptor_node.nunique())
    return MembershipResult(
        receptors,
        rfus,
        edges,
        prepared.unit,
        tuple(identity_cols),
        dict(
            prepared.audit,
            excluded_unassigned_weight=excluded,
            observed_identities=len(receptors),
            assigned_rows=len(original),
            allow_multiple_rfus=allow_multiple_rfus,
            edge_definition="confirmed assignment membership; no sequence or functional distance",
        ),
    )


def prepare_rfu_longitudinal_presence(
    data: Any,
    *,
    sample_col: str,
    sample_order: Sequence[Any],
    identity_cols: Sequence[str] = ("rfu_label",),
    donor_col: str | None = None,
    unit: PlotUnit | None = None,
    min_count: float = 1,
    sample_metadata: pd.DataFrame | None = None,
    missing_identity: str = "exclude",
) -> LongitudinalPresenceResult:
    """Summarize observed identities across explicitly ordered samples/conditions.

    A supplied sample_metadata table declares observed donor/sample combinations,
    including observed zero rows. Without it only combinations present in input
    are observed; absent combinations remain NaN, never imputed as zero. Donor
    identities are kept separate. Omitting donor_col explicitly pools records.
    Missing identity rows can be excluded (audited) or rejected. Denominators use
    included identity-bearing observations in the chosen unit, before detection.
    With two observed samples, flows count detected-union identities, not cells
    or tracked clones. This is descriptive support, not a persistence model.
    """
    if not np.isfinite(min_count) or min_count <= 0:
        raise ValueError("min_count must be finite and positive in the selected observation unit.")
    if missing_identity not in {"exclude", "error"}:
        raise ValueError("missing_identity must be 'exclude' or 'error'.")
    order = _ids(sample_order, context="sample_order").tolist()
    if len(order) < 2:
        raise ValueError("sample_order requires at least two samples/conditions.")
    prepared = _prepared(
        data,
        unit,
        identity_cols if unit == "receptor" and not isinstance(data, ReceptorPlotData) else None,
    )
    if prepared.unit == "receptor" and tuple(identity_cols) != prepared.identity_cols:
        raise ValueError("Longitudinal identities must match receptor-weighted preparation.")
    frame = prepared.rows.copy()
    _require_columns(
        frame,
        [sample_col, *identity_cols, *([donor_col] if donor_col else [])],
        "Longitudinal presence",
    )
    _ids(frame[sample_col], context="Sample labels", unique=False)
    if not frame[sample_col].isin(order).all():
        raise ValueError("sample_order omits observed samples; explicitly subset input first.")
    donor_key = "__plot_donor"
    frame[donor_key] = frame[donor_col] if donor_col else "Pooled observations"
    _ids(frame[donor_key], context="Donor labels", unique=False)
    design = frame[[donor_key, sample_col]].drop_duplicates()
    if sample_metadata is not None:
        _require_columns(
            sample_metadata, [sample_col, *([donor_col] if donor_col else [])], "Sample metadata"
        )
        design = sample_metadata.copy()
        design[donor_key] = design[donor_col] if donor_col else "Pooled observations"
        if (
            design[[donor_key, sample_col]].isna().any().any()
            or design.duplicated([donor_key, sample_col]).any()
        ):
            raise ValueError(
                "Sample metadata requires unique, non-missing donor/sample combinations."
            )
        if not design[sample_col].isin(order).all():
            raise ValueError("Sample metadata contains samples outside sample_order.")
        coverage = frame.merge(
            design[[donor_key, sample_col]], on=[donor_key, sample_col], how="left", indicator=True
        )
        if not coverage._merge.eq("both").all():
            raise ValueError("Input donor/sample combinations are missing from sample_metadata.")
    donors = _ordered(design[donor_key])
    observed = pd.DataFrame(
        False,
        index=_index(donors, name=donor_col or "donor"),
        columns=_index(order, name=sample_col),
    )
    for donor, sample in design[[donor_key, sample_col]].itertuples(index=False, name=None):
        observed.loc[donor, sample] = True
    missing = frame[list(identity_cols)].isna().any(axis=1) | frame[list(identity_cols)].eq("").any(
        axis=1
    )
    excluded = (
        frame.loc[missing]
        .groupby([donor_key, sample_col], observed=True)[_WEIGHT]
        .sum()
        .reset_index()
    )
    if missing.any() and missing_identity == "error":
        raise ValueError("Longitudinal identity fields contain missing values.")
    frame = frame.loc[~missing].copy()
    if frame.empty:
        raise ValueError("No identity-bearing observations remain.")
    frame[_IDENTITY] = _identity(frame, identity_cols)
    frame = _weighted(prepared, frame, [donor_key, sample_col], list(identity_cols))
    den = pd.DataFrame(np.nan, index=observed.index, columns=observed.columns)
    matrices, flows = {}, []
    for donor in donors:
        z = frame.loc[frame[donor_key].eq(donor)]
        identities = sorted(
            z[_IDENTITY].unique(), key=lambda t: tuple((type(x).__name__, str(x)) for x in t)
        )
        matrix = pd.DataFrame(
            np.nan, index=_index(identities, name="identity"), columns=observed.columns
        )
        for sample in order:
            if observed.loc[donor, sample]:
                counts = (
                    z.loc[z[sample_col].eq(sample)].groupby(_IDENTITY, sort=False)[_WEIGHT].sum()
                )
                matrix[sample] = counts.reindex(matrix.index, fill_value=0).to_numpy()
                den.loc[donor, sample] = counts.sum()
        matrices[donor] = matrix
        if len(order) == 2 and observed.loc[donor].all():
            left, right = (matrix[s].ge(min_count) for s in order)
            values = [
                int((left & ~right).sum()),
                int((left & right).sum()),
                int((~left & right).sum()),
            ]
            for status, value in zip(["First only", "Both", "Last only"], values, strict=True):
                flows.append(
                    dict(donor=donor, status=status, identities=value, denominator=sum(values))
                )
    return LongitudinalPresenceResult(
        matrices,
        den,
        observed,
        pd.DataFrame(flows, columns=["donor", "status", "identities", "denominator"]),
        tuple(order),
        tuple(identity_cols),
        prepared.unit,
        float(min_count),
        dict(
            prepared.audit,
            sample_col=sample_col,
            donor_col=donor_col,
            missing_identity=missing_identity,
            excluded_identity_weight=excluded.rename(
                columns={donor_key: "donor", _WEIGHT: "weight"}
            ).to_dict("records"),
            denominator="included identity-bearing weight per donor/sample; paired alluvials use detected identity union",
        ),
    )
