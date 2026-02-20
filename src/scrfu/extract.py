from __future__ import annotations

from typing import Any, Optional

import pandas as pd


def extract_trb_features(
    adata: Any,
    airr_key: str = "airr",
    chain: str = "TRB",
    *,
    prefer_productive: bool = True,
) -> pd.DataFrame:
    """
    Extract per-cell TRB CDR3aa and TRBV.

    Expected (common scirpy pattern):
      adata.obsm[airr_key] is an AIRR-like table (often awkward array or pandas-like).
    We keep this function conservative and support:
      - pandas.DataFrame directly
      - objects convertible to pandas via pd.DataFrame(obj)

    Required columns (case-insensitive handling applied):
      - cell_id (or "cell" / "barcode")
      - chain (or locus)
      - cdr3aa (or junction_aa)
      - v_call (or v_gene)

    Returns a DataFrame with:
      - cell_id
      - cdr3aa
      - trbv
    """
    if not hasattr(adata, "obsm") or airr_key not in adata.obsm:
        raise KeyError(f"AnnData missing obsm['{airr_key}']; cannot extract AIRR features.")

    airr_obj = adata.obsm[airr_key]
    if isinstance(airr_obj, pd.DataFrame):
        df = airr_obj.copy()
    else:
        # best-effort conversion (works for dict-like / record arrays in many cases)
        df = pd.DataFrame(airr_obj)

    # Normalize columns
    cols = {c.lower(): c for c in df.columns}
    def pick(*cands: str) -> Optional[str]:
        for c in cands:
            if c in cols:
                return cols[c]
        return None

    cell_col = pick("cell_id", "cell", "barcode", "cellid")
    chain_col = pick("chain", "locus")
    cdr3_col = pick("cdr3aa", "junction_aa", "cdr3_aa")
    v_col = pick("v_call", "v_gene", "trbv", "v")

    missing = [name for name, col in [
        ("cell_id", cell_col),
        ("chain/locus", chain_col),
        ("cdr3aa/junction_aa", cdr3_col),
        ("v_call/v_gene", v_col),
    ] if col is None]
    if missing:
        raise ValueError(f"AIRR table missing columns: {missing}. Columns seen: {list(df.columns)}")

    # Filter chain
    chain_upper = str(chain).upper()
    chain_series = df[chain_col].astype(str).str.upper()
    df = df.loc[chain_series == chain_upper].copy()

    if df.empty:
        # Return empty but with expected schema
        return pd.DataFrame({"cell_id": [], "cdr3aa": [], "trbv": []})

    # Prefer productive contigs if available
    if prefer_productive:
        prod_col = pick("productive")
        if prod_col is not None:
            # productive may be bool or str
            prod = df[prod_col]
            if prod.dtype != bool:
                prod = prod.astype(str).str.lower().isin(["true", "t", "1", "yes"])
            if prod.any():
                df = df.loc[prod].copy()

    # If there are multiple contigs per cell, choose a stable “top” contig.
    # If a "umis"/"reads" column exists, pick max; else first occurrence.
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

    # Align to adata.obs_names if present: keep only known cells
    if hasattr(adata, "obs_names"):
        obs_names = set(map(str, list(getattr(adata, "obs_names"))))
        out = out[out["cell_id"].isin(obs_names)].copy()

    return out