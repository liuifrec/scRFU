from __future__ import annotations

from typing import Any

import pandas as pd

_COLUMN_ALIASES: dict[str, tuple[str, ...]] = {
    "cell": ("cell_id", "cell", "barcode", "cellid"),
    "chain": ("chain", "locus"),
    "cdr3": ("cdr3aa", "junction_aa", "cdr3_aa", "cdr3"),
    "v": ("v_call", "v_gene", "trbv", "v"),
    "productive": ("productive",),
}


def _n_obs_cells(adata: Any) -> int:
    if hasattr(adata, "obs_names"):
        return len(adata.obs_names)
    if hasattr(adata, "obs"):
        return len(adata.obs)
    return 0


def _obs_names(adata: Any) -> set[str]:
    if hasattr(adata, "obs_names"):
        return set(map(str, list(adata.obs_names)))
    if hasattr(adata, "obs"):
        return set(map(str, list(adata.obs.index)))
    return set()


def _to_dataframe(obj: Any) -> pd.DataFrame:
    return obj.copy() if isinstance(obj, pd.DataFrame) else pd.DataFrame(obj)


def _pick_column(df: pd.DataFrame, aliases: tuple[str, ...]) -> str | None:
    normalized = {str(col).lower(): col for col in df.columns}
    for alias in aliases:
        if alias in normalized:
            return normalized[alias]
    return None


def _truthy_productive(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(False)
    return values.astype(str).str.lower().isin(["true", "t", "1", "yes"])


def _row(
    *,
    airr_key: str,
    chain: str,
    n_obs_cells: int,
    status: str,
    n_airr_rows: int = 0,
    n_airr_cells: int = 0,
    n_overlap_cells: int = 0,
    has_cell_id_col: bool = False,
    has_chain_col: bool = False,
    has_cdr3_col: bool = False,
    has_v_col: bool = False,
    has_productive_col: bool = False,
    n_chain_rows: int = 0,
    n_productive_chain_rows: int = 0,
) -> pd.DataFrame:
    barcode_overlap_rate = n_overlap_cells / n_airr_cells if n_airr_cells else 0.0
    return pd.DataFrame(
        [
            {
                "airr_key": airr_key,
                "chain": chain,
                "n_obs_cells": n_obs_cells,
                "n_airr_rows": n_airr_rows,
                "n_airr_cells": n_airr_cells,
                "n_overlap_cells": n_overlap_cells,
                "barcode_overlap_rate": barcode_overlap_rate,
                "has_cell_id_col": has_cell_id_col,
                "has_chain_col": has_chain_col,
                "has_cdr3_col": has_cdr3_col,
                "has_v_col": has_v_col,
                "has_productive_col": has_productive_col,
                "n_chain_rows": n_chain_rows,
                "n_productive_chain_rows": n_productive_chain_rows,
                "status": status,
            }
        ]
    )


def validate_airr(
    adata: Any,
    airr_key: str = "airr",
    chain: str = "TRB",
    *,
    strict: bool = False,
) -> pd.DataFrame:
    """
    Validate the AIRR-like table used by scRFU without modifying ``adata``.

    The check recognizes common column aliases for cell barcode, chain/locus,
    CDR3 amino-acid sequence, V gene, and productive status. Non-strict mode
    returns a one-row QC table for missing or incomplete inputs. Strict mode
    raises ``ValueError`` when the AIRR key or required columns are missing.
    """
    n_obs = _n_obs_cells(adata)
    if not hasattr(adata, "obsm") or airr_key not in adata.obsm:
        if strict:
            raise ValueError(f"AnnData missing obsm['{airr_key}']")
        return _row(
            airr_key=airr_key,
            chain=chain,
            n_obs_cells=n_obs,
            status="missing_airr_key",
        )

    try:
        df = _to_dataframe(adata.obsm[airr_key])
    except Exception as exc:
        if strict:
            raise ValueError(f"Could not convert obsm['{airr_key}'] to a table") from exc
        return _row(airr_key=airr_key, chain=chain, n_obs_cells=n_obs, status="invalid_airr")

    cell_col = _pick_column(df, _COLUMN_ALIASES["cell"])
    chain_col = _pick_column(df, _COLUMN_ALIASES["chain"])
    cdr3_col = _pick_column(df, _COLUMN_ALIASES["cdr3"])
    v_col = _pick_column(df, _COLUMN_ALIASES["v"])
    productive_col = _pick_column(df, _COLUMN_ALIASES["productive"])

    has_cell = cell_col is not None
    has_chain = chain_col is not None
    has_cdr3 = cdr3_col is not None
    has_v = v_col is not None
    has_productive = productive_col is not None

    missing = [
        name
        for name, present in [
            ("cell_id", has_cell),
            ("chain/locus", has_chain),
            ("cdr3aa/junction_aa", has_cdr3),
            ("v_call/v_gene", has_v),
        ]
        if not present
    ]
    if missing:
        n_airr_cells = 0
        n_overlap_cells = 0
        if has_cell:
            airr_cells = set(df[cell_col].dropna().astype(str))
            n_airr_cells = len(airr_cells)
            n_overlap_cells = len(airr_cells.intersection(_obs_names(adata)))

        n_chain_rows = 0
        n_productive_chain_rows = 0
        if has_chain:
            chain_mask = df[chain_col].astype(str).str.upper() == str(chain).upper()
            n_chain_rows = int(chain_mask.sum())
            if has_productive:
                productive = _truthy_productive(df[productive_col])
                n_productive_chain_rows = int((chain_mask & productive).sum())

        if strict:
            raise ValueError(f"AIRR table missing required columns: {missing}")
        return _row(
            airr_key=airr_key,
            chain=chain,
            n_obs_cells=n_obs,
            n_airr_rows=len(df),
            n_airr_cells=n_airr_cells,
            n_overlap_cells=n_overlap_cells,
            has_cell_id_col=has_cell,
            has_chain_col=has_chain,
            has_cdr3_col=has_cdr3,
            has_v_col=has_v,
            has_productive_col=has_productive,
            n_chain_rows=n_chain_rows,
            n_productive_chain_rows=n_productive_chain_rows,
            status="missing_required_columns",
        )

    cell_values = df[cell_col].dropna().astype(str)
    airr_cells = set(cell_values)
    obs_names = _obs_names(adata)
    overlap = airr_cells.intersection(obs_names)

    chain_upper = str(chain).upper()
    chain_mask = df[chain_col].astype(str).str.upper() == chain_upper
    n_chain_rows = int(chain_mask.sum())

    n_productive_chain_rows = 0
    if has_productive:
        productive = _truthy_productive(df[productive_col])
        n_productive_chain_rows = int((chain_mask & productive).sum())

    status = "ok"
    if n_chain_rows == 0:
        status = "no_chain_rows"
    elif airr_cells and not overlap:
        status = "barcode_mismatch"

    return _row(
        airr_key=airr_key,
        chain=chain,
        n_obs_cells=n_obs,
        n_airr_rows=len(df),
        n_airr_cells=len(airr_cells),
        n_overlap_cells=len(overlap),
        has_cell_id_col=has_cell,
        has_chain_col=has_chain,
        has_cdr3_col=has_cdr3,
        has_v_col=has_v,
        has_productive_col=has_productive,
        n_chain_rows=n_chain_rows,
        n_productive_chain_rows=n_productive_chain_rows,
        status=status,
    )


__all__ = ["validate_airr"]
