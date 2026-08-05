from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import pandas as pd

_WELLS_TCR_IR_KEY = "TCR_IR"
_WELLS_COLUMN_ALIASES: dict[str, tuple[str, ...]] = {
    "cell": ("cell_id", "cell", "barcode", "cellid"),
    "chain": ("tcr-ir_vdj_1_locus", "ir_vdj_1_locus"),
    "cdr3": (
        "tcr-ir_vdj_1_junction_aa",
        "ir_vdj_1_junction_aa",
        "tcr-ir_vdj_1_cdr3",
        "ir_vdj_1_cdr3",
    ),
    "v": ("tcr-ir_vdj_1_v_call", "ir_vdj_1_v_call"),
    "productive": ("tcr-ir_vdj_1_productive", "ir_vdj_1_productive"),
}


def _pick_column(df: pd.DataFrame, aliases: tuple[str, ...]) -> Any | None:
    normalized = {str(column).lower(): column for column in df.columns}
    for alias in aliases:
        if alias in normalized:
            return normalized[alias]
    return None


def _truthy(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(False)
    return values.astype(str).str.lower().isin(["true", "t", "1", "yes"])


def _present_text(values: pd.Series) -> pd.Series:
    normalized = values.astype("string").str.strip().str.lower()
    return values.notna() & ~normalized.isin(["", "nan", "none", "na", "<na>"])


def _obs_names(adata: Any) -> list[str]:
    if not hasattr(adata, "obs_names"):
        raise ValueError("Wells TCR_IR extraction requires adata.obs_names for row alignment.")
    return list(map(str, list(adata.obs_names)))


def _wells_tcr_ir_table(adata: Any, key: str) -> tuple[pd.DataFrame, str]:
    """Return the flattened Wells TCR table and its AnnData container."""
    if hasattr(adata, "uns") and key in adata.uns:
        obj = adata.uns[key]
        container = "uns"
    elif hasattr(adata, "obsm") and key in adata.obsm:
        obj = adata.obsm[key]
        container = "obsm"
    else:
        raise KeyError(
            f"AnnData missing uns['{key}'] and obsm['{key}']; cannot extract Wells TCR_IR features."
        )

    if isinstance(obj, pd.DataFrame):
        return obj.copy(), container
    if isinstance(obj, Mapping):
        return pd.DataFrame(dict(obj)), container
    try:
        return pd.DataFrame(obj), container
    except Exception as exc:
        raise ValueError(f"Could not convert {container}['{key}'] to a table.") from exc


def _wells_cell_ids(adata: Any, df: pd.DataFrame, cell_col: Any | None) -> pd.Series:
    """Resolve cell IDs without silently guessing a non-positional index."""
    if cell_col is not None:
        cells = df[cell_col]
    else:
        obs_names = _obs_names(adata)
        if len(df) != len(obs_names):
            raise ValueError(
                "Wells TCR_IR table has no cell column and cannot be aligned to adata.obs_names: "
                f"{len(df)} table rows versus {len(obs_names)} observations."
            )

        index_as_str = list(map(str, df.index))
        default_index = isinstance(df.index, pd.RangeIndex) and df.index.equals(
            pd.RangeIndex(len(df))
        )
        if not default_index and index_as_str != obs_names:
            raise ValueError(
                "Wells TCR_IR table has no cell column and its index does not exactly match "
                "adata.obs_names; refusing ambiguous row alignment."
            )
        cells = pd.Series(obs_names, index=df.index)

    cell_strings = cells.astype("string")
    valid_cells = cell_strings.dropna()
    if valid_cells.duplicated().any():
        duplicates = sorted(set(valid_cells[valid_cells.duplicated(keep=False)].astype(str)))
        raise ValueError(f"Wells TCR_IR table contains duplicate cell identifiers: {duplicates}")
    return cell_strings


def extract_wells_tcr_ir_features(
    adata: Any,
    key: str = _WELLS_TCR_IR_KEY,
    chain: str = "TRB",
    *,
    prefer_productive: bool = True,
) -> pd.DataFrame:
    """Adapt the Wells atlas flattened ``TCR_IR`` table to scRFU features.

    The Wells atlas processing workflow flattened Scirpy's primary VDJ fields
    into columns such as ``TCR-IR_VDJ_1_junction_aa`` and stored the resulting
    table in ``obsm["TCR_IR"]``. Some distributed objects place the same table
    in ``uns["TCR_IR"]``. When both are present, ``uns`` takes precedence.

    Cell identifiers must be explicit or exactly row-aligned to
    ``adata.obs_names``. The adapter does not mutate ``adata``.
    """
    df, container = _wells_tcr_ir_table(adata, key)

    cell_col = _pick_column(df, _WELLS_COLUMN_ALIASES["cell"])
    chain_col = _pick_column(df, _WELLS_COLUMN_ALIASES["chain"])
    cdr3_col = _pick_column(df, _WELLS_COLUMN_ALIASES["cdr3"])
    v_col = _pick_column(df, _WELLS_COLUMN_ALIASES["v"])
    productive_col = _pick_column(df, _WELLS_COLUMN_ALIASES["productive"])

    missing = [
        label
        for label, column in [("VDJ_1 CDR3 amino acid", cdr3_col), ("VDJ_1 V call", v_col)]
        if column is None
    ]
    if missing:
        raise ValueError(
            f"Wells {container}['{key}'] table missing columns: {missing}. "
            f"Columns seen: {list(df.columns)}"
        )

    cell_ids = _wells_cell_ids(adata, df, cell_col)
    keep = pd.Series(True, index=df.index)

    if chain_col is not None:
        keep &= df[chain_col].astype("string").str.upper().eq(str(chain).upper()).fillna(False)

    if prefer_productive and productive_col is not None:
        productive = _truthy(df[productive_col])
        if productive.any():
            keep &= productive

    cdr3 = df[cdr3_col].astype("string").str.strip()
    trbv = df[v_col].astype("string").str.strip()
    keep &= _present_text(cdr3) & _present_text(trbv)

    obs_names = set(_obs_names(adata))
    keep &= cell_ids.isin(obs_names)

    out = pd.DataFrame(
        {
            "cell_id": cell_ids.loc[keep].astype(str).to_numpy(),
            "cdr3aa": cdr3.loc[keep].astype(str).to_numpy(),
            "trbv": trbv.loc[keep].astype(str).to_numpy(),
        }
    )
    return out.reset_index(drop=True)


def extract_trb_features(
    adata: Any,
    airr_key: str = "airr",
    chain: str = "TRB",
    *,
    prefer_productive: bool = True,
) -> pd.DataFrame:
    """
    Extract per-cell TRB CDR3aa and TRBV.

    The usual input is an AIRR-like table in ``adata.obsm[airr_key]``. Passing
    ``airr_key="TCR_IR"`` also recognizes the Wells atlas flattened table in
    ``adata.uns["TCR_IR"]`` or ``adata.obsm["TCR_IR"]``.

    Returns a DataFrame with ``cell_id``, ``cdr3aa``, and ``trbv``.
    """
    if airr_key == _WELLS_TCR_IR_KEY and (
        (hasattr(adata, "uns") and airr_key in adata.uns)
        or (hasattr(adata, "obsm") and airr_key in adata.obsm)
    ):
        return extract_wells_tcr_ir_features(
            adata,
            key=airr_key,
            chain=chain,
            prefer_productive=prefer_productive,
        )

    if not hasattr(adata, "obsm") or airr_key not in adata.obsm:
        raise KeyError(f"AnnData missing obsm['{airr_key}']; cannot extract AIRR features.")

    airr_obj = adata.obsm[airr_key]
    df = airr_obj.copy() if isinstance(airr_obj, pd.DataFrame) else pd.DataFrame(airr_obj)

    cols = {str(c).lower(): c for c in df.columns}

    def pick(*cands: str) -> Any | None:
        for candidate in cands:
            if candidate in cols:
                return cols[candidate]
        return None

    cell_col = pick("cell_id", "cell", "barcode", "cellid")
    chain_col = pick("chain", "locus")
    cdr3_col = pick("cdr3aa", "junction_aa", "cdr3_aa", "cdr3")
    v_col = pick("v_call", "v_gene", "trbv", "v")

    missing = [
        name
        for name, col in [
            ("cell_id", cell_col),
            ("chain/locus", chain_col),
            ("cdr3aa/junction_aa", cdr3_col),
            ("v_call/v_gene", v_col),
        ]
        if col is None
    ]
    if missing:
        raise ValueError(f"AIRR table missing columns: {missing}. Columns seen: {list(df.columns)}")

    chain_upper = str(chain).upper()
    chain_series = df[chain_col].astype(str).str.upper()
    df = df.loc[chain_series == chain_upper].copy()

    if df.empty:
        return pd.DataFrame({"cell_id": [], "cdr3aa": [], "trbv": []})

    if prefer_productive:
        prod_col = pick("productive")
        if prod_col is not None:
            prod = _truthy(df[prod_col])
            if prod.any():
                df = df.loc[prod].copy()

    umi_col = pick("umis", "umi", "reads", "read_count", "junction_reads")
    if umi_col is not None:
        df["_rank"] = pd.to_numeric(df[umi_col], errors="coerce").fillna(0)
        df = df.sort_values(["_rank"], ascending=False).drop(columns=["_rank"])
    df = df.drop_duplicates(subset=[cell_col], keep="first")

    out = pd.DataFrame(
        {
            "cell_id": df[cell_col].astype(str).values,
            "cdr3aa": df[cdr3_col].astype(str).values,
            "trbv": df[v_col].astype(str).values,
        }
    )

    if hasattr(adata, "obs_names"):
        obs_names = set(map(str, list(adata.obs_names)))
        out = out[out["cell_id"].isin(obs_names)].copy()

    return out


__all__ = ["extract_trb_features", "extract_wells_tcr_ir_features"]
