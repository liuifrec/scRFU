from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import pandas as pd

RECEPTOR_SCHEMA_VERSION = "1.0"
REQUIRED_RECEPTOR_COLUMNS = (
    "input_row_id",
    "cell_id",
    "chain",
    "cdr3aa",
    "v_call",
    "productive",
    "source_adapter",
    "source_row_id",
)
OPTIONAL_RECEPTOR_COLUMNS = (
    "duplicate_count",
    "consensus_count",
    "source_slot",
    "j_call",
    "d_call",
    "c_call",
    "junction",
    "junction_aa",
    "umi_count",
    "read_count",
    "clonotype_id",
    "receptor_type",
    "light_chain",
)
RECOGNIZED_CHAINS = frozenset({"TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL"})

_CHAIN_ALIASES = {
    "A": "TRA",
    "ALPHA": "TRA",
    "TCRA": "TRA",
    "B": "TRB",
    "BETA": "TRB",
    "TCRB": "TRB",
    "G": "TRG",
    "GAMMA": "TRG",
    "TCRG": "TRG",
    "D": "TRD",
    "DELTA": "TRD",
    "TCRD": "TRD",
    "HEAVY": "IGH",
    "KAPPA": "IGK",
    "LAMBDA": "IGL",
}
_TRUE = frozenset({"true", "t", "1", "yes", "y", "productive"})
_FALSE = frozenset({"false", "f", "0", "no", "n", "nonproductive", "unproductive"})
_MISSING_TEXT = frozenset({"", "nan", "none", "na", "<na>"})


def normalize_chain(value: Any) -> Any:
    if pd.isna(value):
        return pd.NA
    normalized = str(value).strip().upper().replace("-", "").replace("_", "")
    if normalized in _MISSING_TEXT:
        return pd.NA
    return _CHAIN_ALIASES.get(normalized, normalized)


def normalize_productive(value: Any) -> Any:
    if pd.isna(value):
        return pd.NA
    if isinstance(value, bool):
        return value
    if isinstance(value, int) and value in {0, 1}:
        return bool(value)
    normalized = str(value).strip().lower()
    if normalized in _TRUE:
        return True
    if normalized in _FALSE:
        return False
    if normalized in _MISSING_TEXT:
        return pd.NA
    return pd.NA


def normalize_text(value: Any) -> Any:
    if pd.isna(value):
        return pd.NA
    normalized = str(value).strip()
    if normalized.lower() in _MISSING_TEXT:
        return pd.NA
    return normalized


def canonicalize_receptor_table(
    table: pd.DataFrame,
    *,
    source_adapter: str | None = None,
    source_row_ids: pd.Series | None = None,
) -> pd.DataFrame:
    """Return a normalized canonical receptor table without changing source order."""
    if not isinstance(table, pd.DataFrame):
        raise TypeError("receptor table must be a pandas DataFrame.")
    out = table.copy().reset_index(drop=True)
    if "input_row_id" not in out:
        out.insert(0, "input_row_id", [f"row_{index:08d}" for index in range(len(out))])
    if "source_row_id" not in out:
        values = source_row_ids if source_row_ids is not None else table.index.to_series(index=None)
        out["source_row_id"] = pd.Series(values, dtype="string").reset_index(drop=True)
    if "source_adapter" not in out and source_adapter is not None:
        out["source_adapter"] = source_adapter
    if "productive" not in out:
        out["productive"] = pd.Series(pd.NA, index=out.index, dtype="boolean")

    for column in ("cell_id", "source_adapter", "source_row_id"):
        if column in out:
            out[column] = out[column].map(normalize_text).astype("string")
    if "chain" in out:
        out["chain"] = out["chain"].map(normalize_chain).astype("string")
    for column in ("cdr3aa", "v_call"):
        if column in out:
            out[column] = out[column].map(normalize_text).astype("string")
    out["productive"] = pd.array(out["productive"].map(normalize_productive), dtype="boolean")

    first = [column for column in REQUIRED_RECEPTOR_COLUMNS if column in out]
    remainder = [column for column in out.columns if column not in first]
    return out.loc[:, [*first, *remainder]]


def validate_receptor_table(table: pd.DataFrame, strict: bool = True) -> dict[str, Any]:
    """Validate the dataset-independent receptor-row contract and return actionable QC."""
    if not isinstance(table, pd.DataFrame):
        raise TypeError("receptor table must be a pandas DataFrame.")
    missing_columns = [column for column in REQUIRED_RECEPTOR_COLUMNS if column not in table]
    errors: list[str] = []
    if missing_columns:
        errors.append(f"missing required columns: {missing_columns}")
        normalized = table.copy()
    else:
        normalized = canonicalize_receptor_table(table)

    def missing_count(column: str) -> int:
        return int(normalized[column].isna().sum()) if column in normalized else len(normalized)

    duplicated_input = (
        normalized.loc[normalized["input_row_id"].duplicated(keep=False), "input_row_id"]
        .dropna()
        .astype(str)
        .tolist()
        if "input_row_id" in normalized
        else []
    )
    duplicated_source = (
        normalized.loc[normalized["source_row_id"].duplicated(keep=False), "source_row_id"]
        .dropna()
        .astype(str)
        .tolist()
        if "source_row_id" in normalized
        else []
    )
    if missing_count("input_row_id"):
        errors.append("input_row_id contains missing values")
    if duplicated_input:
        errors.append(f"input_row_id contains duplicates: {sorted(set(duplicated_input))[:10]}")
    if missing_count("cell_id"):
        errors.append("cell_id contains missing values")
    if missing_count("cdr3aa"):
        errors.append("cdr3aa contains missing values")
    if missing_count("source_adapter"):
        errors.append("source_adapter contains missing values")
    if missing_count("source_row_id"):
        errors.append("source_row_id contains missing values")

    chain_values = (
        sorted(normalized["chain"].dropna().astype(str).unique().tolist())
        if "chain" in normalized
        else []
    )
    unrecognized = sorted(set(chain_values).difference(RECOGNIZED_CHAINS))
    if strict and unrecognized:
        errors.append(f"unrecognized chain values: {unrecognized}")

    report: dict[str, Any] = {
        "schema_version": RECEPTOR_SCHEMA_VERSION,
        "status": "ok" if not errors else "error",
        "row_count": len(table),
        "unique_cell_count": int(normalized["cell_id"].nunique(dropna=True))
        if "cell_id" in normalized
        else 0,
        "unique_input_row_id_count": int(normalized["input_row_id"].nunique(dropna=True))
        if "input_row_id" in normalized
        else 0,
        "missing_cell_id_count": missing_count("cell_id"),
        "missing_cdr3aa_count": missing_count("cdr3aa"),
        "recognized_chains": sorted(set(chain_values).intersection(RECOGNIZED_CHAINS)),
        "unrecognized_chains": unrecognized,
        "productive_value_counts": {
            str(key): int(value)
            for key, value in normalized.get("productive", pd.Series(dtype="boolean"))
            .value_counts(dropna=False)
            .items()
        },
        "duplicated_input_row_ids": sorted(set(duplicated_input)),
        "duplicated_source_row_ids": sorted(set(duplicated_source)),
        "errors": errors,
    }
    if strict and errors:
        raise ValueError("Invalid canonical receptor table: " + "; ".join(errors))
    return report


def receptor_schema() -> Mapping[str, Any]:
    return {
        "version": RECEPTOR_SCHEMA_VERSION,
        "required_columns": list(REQUIRED_RECEPTOR_COLUMNS),
        "optional_columns": list(OPTIONAL_RECEPTOR_COLUMNS),
        "recognized_chains": sorted(RECOGNIZED_CHAINS),
    }


__all__ = [
    "OPTIONAL_RECEPTOR_COLUMNS",
    "RECEPTOR_SCHEMA_VERSION",
    "RECOGNIZED_CHAINS",
    "REQUIRED_RECEPTOR_COLUMNS",
    "canonicalize_receptor_table",
    "normalize_chain",
    "normalize_productive",
    "receptor_schema",
    "validate_receptor_table",
]
