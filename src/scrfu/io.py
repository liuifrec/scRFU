from __future__ import annotations

import hashlib
import json
import platform
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from ._version import __version__
from .pp import RECEPTOR_SCHEMA_VERSION, validate_receptor_table
from .summary import aggregate_rfu

PathLike = str | Path
RECEPTOR_CACHE_SCHEMA_VERSION = 2


class UnsupportedH5ADLayout(ValueError):
    """Raised when selective reading encounters an unsupported H5AD encoding."""


@dataclass(frozen=True)
class H5ADDataFrameInfo:
    location: str
    key: str
    row_count: int
    index_name: str
    columns: tuple[str, ...]


@dataclass(frozen=True)
class ReceptorCacheData:
    receptors: pd.DataFrame
    cell_metadata: pd.DataFrame
    manifest: dict[str, Any]
    validation: dict[str, Any]


def read_h5ad(path: PathLike) -> ad.AnnData:
    return ad.read_h5ad(str(path))


def write_h5ad(adata: ad.AnnData, path: PathLike) -> None:
    adata.write_h5ad(str(path))


def _text(value: Any) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def _encoding_type(element: Any) -> str:
    return _text(element.attrs.get("encoding-type", ""))


def _encoding_version(element: Any) -> str:
    return _text(element.attrs.get("encoding-version", ""))


def _string_list(values: Any) -> list[str]:
    if isinstance(values, (str, bytes)):
        return [_text(values)]
    return [_text(value) for value in values]


def _validate_dataframe(group: Any, logical_name: str) -> tuple[str, list[str], int]:
    import h5py

    if not isinstance(group, h5py.Group) or _encoding_type(group) != "dataframe":
        seen = _encoding_type(group) or "missing"
        raise UnsupportedH5ADLayout(
            f"{logical_name} has unsupported encoding-type {seen!r}; expected 'dataframe'."
        )
    version = _encoding_version(group)
    if version != "0.2.0":
        raise UnsupportedH5ADLayout(
            f"{logical_name} has unsupported dataframe encoding-version {version or 'missing'!r}."
        )
    index_key = _text(group.attrs.get("_index", ""))
    if not index_key or index_key not in group:
        raise UnsupportedH5ADLayout(f"{logical_name} does not declare a readable '_index'.")
    index = group[index_key]
    if not isinstance(index, h5py.Dataset) or index.ndim != 1:
        raise UnsupportedH5ADLayout(f"{logical_name} index must be a one-dimensional dataset.")
    if _encoding_type(index) not in {"array", "string-array"}:
        raise UnsupportedH5ADLayout(
            f"{logical_name} index has unsupported encoding-type "
            f"{_encoding_type(index) or 'missing'!r}."
        )
    if "column-order" not in group.attrs:
        raise UnsupportedH5ADLayout(f"{logical_name} is missing 'column-order'.")
    columns = _string_list(group.attrs["column-order"])
    missing = [column for column in columns if column not in group]
    if missing:
        raise UnsupportedH5ADLayout(f"{logical_name} declares missing columns: {missing}.")
    row_count = int(index.shape[0])
    for column in columns:
        element = group[column]
        if isinstance(element, h5py.Dataset):
            column_rows = element.shape[0] if element.ndim else -1
        elif isinstance(element, h5py.Group) and "codes" in element:
            column_rows = element["codes"].shape[0]
        elif isinstance(element, h5py.Group) and "values" in element:
            column_rows = element["values"].shape[0]
        else:
            column_rows = row_count
        if column_rows != row_count:
            raise UnsupportedH5ADLayout(
                f"{logical_name}[{column!r}] has {column_rows} rows; expected {row_count}."
            )
    return index_key, columns, row_count


def _decode(values: np.ndarray) -> np.ndarray:
    if values.dtype.kind in {"S", "O", "U"}:
        return np.asarray(
            [value.decode("utf-8") if isinstance(value, bytes) else value for value in values],
            dtype=object,
        )
    return values


def _take(dataset: Any, positions: np.ndarray) -> np.ndarray:
    if len(positions) == 0:
        return np.asarray(dataset[:0])
    order = np.argsort(positions, kind="stable")
    sorted_positions = positions[order]
    if len(np.unique(sorted_positions)) != len(sorted_positions):
        raise ValueError("Selected H5AD positions must be unique.")
    values = np.asarray(dataset[sorted_positions])
    inverse = np.empty_like(order)
    inverse[order] = np.arange(len(order))
    return values[inverse]


def _read_column(element: Any, positions: np.ndarray, logical_name: str) -> Any:
    import h5py

    if isinstance(element, h5py.Dataset):
        encoding = _encoding_type(element)
        if encoding not in {"array", "string-array"} or element.ndim != 1:
            raise UnsupportedH5ADLayout(
                f"{logical_name} has unsupported dataset encoding {encoding or 'missing'!r}."
            )
        return _decode(_take(element, positions))
    if not isinstance(element, h5py.Group):
        raise UnsupportedH5ADLayout(f"{logical_name} is not a supported HDF5 element.")
    encoding = _encoding_type(element)
    if encoding == "categorical":
        if _encoding_version(element) != "0.2.0" or not {"codes", "categories"}.issubset(element):
            raise UnsupportedH5ADLayout(f"{logical_name} has an invalid categorical encoding.")
        codes = _take(element["codes"], positions).astype(np.int64, copy=False)
        categories = _decode(np.asarray(element["categories"][:]))
        if np.any(codes < -1) or np.any(codes >= len(categories)):
            raise UnsupportedH5ADLayout(f"{logical_name} categorical codes are out of range.")
        return pd.Categorical.from_codes(
            codes, categories=categories, ordered=bool(element.attrs.get("ordered", False))
        )
    if encoding in {"nullable-integer", "nullable-boolean", "nullable-string-array"}:
        if not {"values", "mask"}.issubset(element):
            raise UnsupportedH5ADLayout(f"{logical_name} has an invalid {encoding} encoding.")
        values = _decode(_take(element["values"], positions))
        mask = _take(element["mask"], positions).astype(bool, copy=False)
        series = pd.Series(values)
        series.loc[mask] = pd.NA
        dtype = {
            "nullable-integer": "Int64",
            "nullable-boolean": "boolean",
            "nullable-string-array": "string",
        }[encoding]
        return pd.array(series, dtype=dtype)
    raise UnsupportedH5ADLayout(
        f"{logical_name} has unsupported encoding-type {encoding or 'missing'!r}."
    )


def _validate_h5ad_root(handle: Any) -> int:
    if _encoding_type(handle) != "anndata":
        raise UnsupportedH5ADLayout(
            f"H5AD root has encoding-type {_encoding_type(handle) or 'missing'!r}; expected 'anndata'."
        )
    if _encoding_version(handle) not in {"0.1.0"}:
        raise UnsupportedH5ADLayout(
            f"H5AD root has unsupported encoding-version "
            f"{_encoding_version(handle) or 'missing'!r}."
        )
    if "obs" not in handle:
        raise UnsupportedH5ADLayout("H5AD root is missing the obs dataframe.")
    _, _, obs_count = _validate_dataframe(handle["obs"], "obs")
    return obs_count


def _positions_for_index(
    group: Any,
    *,
    row_count: int,
    selected_names: Sequence[str] | None,
    max_rows: int | None,
) -> np.ndarray:
    if selected_names is not None and max_rows is not None:
        raise ValueError("Pass selected_names or max_rows, not both.")
    if max_rows is not None:
        if isinstance(max_rows, bool) or not isinstance(max_rows, int) or max_rows <= 0:
            raise ValueError("max_rows must be a positive integer.")
        return np.arange(min(max_rows, row_count), dtype=np.int64)
    if selected_names is None:
        return np.arange(row_count, dtype=np.int64)
    index_key, _, _ = _validate_dataframe(group, "selected dataframe")
    names = [_text(value) for value in _decode(np.asarray(group[index_key][:]))]
    if len(set(names)) != len(names):
        raise UnsupportedH5ADLayout("Selected dataframe index contains duplicates.")
    lookup = {name: position for position, name in enumerate(names)}
    requested = list(map(str, selected_names))
    missing = [name for name in requested if name not in lookup]
    if missing:
        raise KeyError(f"Selected row names are absent from the H5AD: {missing[:10]}")
    return np.asarray([lookup[name] for name in requested], dtype=np.int64)


def _read_dataframe_group(
    group: Any,
    *,
    logical_name: str,
    columns: Sequence[str] | None,
    positions: np.ndarray,
) -> tuple[pd.DataFrame, H5ADDataFrameInfo]:
    index_key, available, row_count = _validate_dataframe(group, logical_name)
    if np.any(positions < 0) or np.any(positions >= row_count):
        raise IndexError(f"Selected positions fall outside {logical_name} ({row_count} rows).")
    selected = available if columns is None else list(columns)
    missing = [column for column in selected if column not in available]
    if missing:
        raise ValueError(
            f"Requested columns are missing from {logical_name}: {missing}. Available: {available}"
        )
    index = pd.Index(
        [_text(value) for value in _decode(_take(group[index_key], positions))], name=None
    )
    data = {
        column: _read_column(group[column], positions, f"{logical_name}[{column!r}]")
        for column in selected
    }
    frame = pd.DataFrame(data, index=index)
    info = H5ADDataFrameInfo(
        location=logical_name.split("[")[0],
        key=logical_name,
        row_count=row_count,
        index_name=index_key,
        columns=tuple(available),
    )
    return frame, info


def read_h5ad_dataframe(
    path: PathLike,
    *,
    location: str = "uns",
    key: str,
    selected_names: Sequence[str] | None = None,
    max_rows: int | None = None,
    return_info: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, H5ADDataFrameInfo]:
    """Selectively read one encoded H5AD dataframe without materializing expression data."""
    import h5py

    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"H5AD not found: {source}")
    if location not in {"uns", "obsm"}:
        raise ValueError("location must be 'uns' or 'obsm'.")
    with h5py.File(source, "r") as handle:
        obs_count = _validate_h5ad_root(handle)
        hdf5_key = f"{location}/{key}"
        if hdf5_key not in handle:
            raise KeyError(f"H5AD is missing {location}[{key!r}].")
        group = handle[hdf5_key]
        _, _, row_count = _validate_dataframe(group, f"{location}[{key!r}]")
        if location == "obsm" and row_count != obs_count:
            raise UnsupportedH5ADLayout(
                f"obsm[{key!r}] has {row_count} rows but obs has {obs_count}."
            )
        if location == "obsm":
            obsm_index_key, _, _ = _validate_dataframe(group, f"obsm[{key!r}]")
            obs_index_key, _, _ = _validate_dataframe(handle["obs"], "obs")
            obsm_index = [_text(value) for value in _decode(np.asarray(group[obsm_index_key][:]))]
            obs_index = [
                _text(value) for value in _decode(np.asarray(handle["obs"][obs_index_key][:]))
            ]
            if obsm_index != obs_index:
                raise UnsupportedH5ADLayout(
                    f"obsm[{key!r}] index does not exactly match the obs index."
                )
        positions = _positions_for_index(
            group, row_count=row_count, selected_names=selected_names, max_rows=max_rows
        )
        frame, info = _read_dataframe_group(
            group, logical_name=f"{location}[{key!r}]", columns=None, positions=positions
        )
    return (frame, info) if return_info else frame


def read_h5ad_obs(
    path: PathLike,
    *,
    columns: Sequence[str] | None = None,
    selected_names: Sequence[str] | None = None,
    max_rows: int | None = None,
    return_info: bool = False,
) -> pd.DataFrame | tuple[pd.DataFrame, H5ADDataFrameInfo]:
    """Read selected observation metadata and preserve exact H5AD observation order."""
    import h5py

    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"H5AD not found: {source}")
    with h5py.File(source, "r") as handle:
        obs_count = _validate_h5ad_root(handle)
        group = handle["obs"]
        positions = _positions_for_index(
            group, row_count=obs_count, selected_names=selected_names, max_rows=max_rows
        )
        frame, info = _read_dataframe_group(
            group, logical_name="obs", columns=columns, positions=positions
        )
    return (frame, info) if return_info else frame


def read_h5ad_shape(path: PathLike) -> tuple[int, int]:
    """Read and validate H5AD dimensions from metadata without reading matrix values."""
    import h5py

    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"H5AD not found: {source}")
    with h5py.File(source, "r") as handle:
        obs_count = _validate_h5ad_root(handle)
        if "X" in handle:
            matrix = handle["X"]
            if isinstance(matrix, h5py.Dataset):
                shape = tuple(map(int, matrix.shape))
            elif isinstance(matrix, h5py.Group) and "shape" in matrix.attrs:
                shape = tuple(map(int, matrix.attrs["shape"]))
            else:
                raise UnsupportedH5ADLayout(
                    "H5AD X has neither a dataset shape nor a sparse-group shape attribute."
                )
        elif "var" in handle:
            _, _, var_count = _validate_dataframe(handle["var"], "var")
            shape = (obs_count, var_count)
        else:
            raise UnsupportedH5ADLayout("H5AD has neither X nor a var dataframe.")
        if len(shape) != 2 or shape[0] != obs_count:
            raise UnsupportedH5ADLayout(
                f"H5AD X shape {shape} is inconsistent with {obs_count} observation rows."
            )
        return shape


def file_sha256(path: PathLike) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_fingerprint(path: PathLike) -> dict[str, Any]:
    source = Path(path).expanduser().resolve()
    stat = source.stat()
    sample_size = 1024 * 1024
    offsets = sorted(
        {0, max(0, stat.st_size // 2 - sample_size // 2), max(0, stat.st_size - sample_size)}
    )
    digest = hashlib.sha256()
    digest.update(f"size:{stat.st_size}\n".encode())
    with source.open("rb") as stream:
        for offset in offsets:
            stream.seek(offset)
            block = stream.read(sample_size)
            digest.update(f"offset:{offset}:length:{len(block)}\n".encode())
            digest.update(block)
    return {
        "path": str(source),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "fingerprint_method": "sha256-size-first-middle-last-1MiB-v1",
        "fingerprint": digest.hexdigest(),
    }


def _metadata_for_cache(metadata: pd.DataFrame) -> pd.DataFrame:
    result = metadata.copy()
    if "cell_id" not in result:
        result.insert(0, "cell_id", result.index.astype(str))
    result["cell_id"] = result["cell_id"].astype("string")
    if result["cell_id"].isna().any() or result["cell_id"].duplicated().any():
        raise ValueError("Cell metadata must contain unique, non-missing cell_id values.")
    return result.reset_index(drop=True)


def write_receptor_cache(
    cache_dir: PathLike,
    receptors: pd.DataFrame,
    cell_metadata: pd.DataFrame,
    *,
    source_adapter: str,
    source_adapter_version: str,
    source_format: str,
    source_path: PathLike | None = None,
    adapter_qc: dict[str, Any] | None = None,
    adapter_configuration: dict[str, Any] | None = None,
    selected_metadata_columns: Sequence[str] = (),
    source_atlas_dimensions: Sequence[int] | None = None,
    force: bool = False,
) -> dict[str, Any]:
    """Write a portable checksummed canonical receptor cache."""
    validate_receptor_table(receptors, strict=True)
    metadata = _metadata_for_cache(cell_metadata)
    if len(metadata):
        unknown = sorted(
            set(receptors["cell_id"].astype(str)).difference(metadata["cell_id"].astype(str))
        )
        if unknown:
            raise ValueError(f"Receptor cells are absent from cell metadata: {unknown[:10]}")
    output = Path(cache_dir).expanduser().resolve()
    existing = [
        output / name
        for name in ("receptors.tsv.gz", "obs_metadata.tsv.gz", "preparation_manifest.json")
    ]
    if not force and any(path.exists() for path in existing):
        raise FileExistsError(
            f"Receptor cache already exists at {output}; pass force=True to replace it."
        )
    output.mkdir(parents=True, exist_ok=True)
    receptor_path = output / "receptors.tsv.gz"
    metadata_path = output / "obs_metadata.tsv.gz"
    receptors.to_csv(receptor_path, sep="\t", index=False, compression="gzip")
    metadata.to_csv(metadata_path, sep="\t", index=False, compression="gzip")
    source = source_fingerprint(source_path) if source_path is not None else None
    manifest = {
        "cache_schema_version": RECEPTOR_CACHE_SCHEMA_VERSION,
        "receptor_schema_version": RECEPTOR_SCHEMA_VERSION,
        "source_adapter": source_adapter,
        "source_adapter_version": source_adapter_version,
        "source_format": source_format,
        "source": source,
        "source_atlas_dimensions": list(source_atlas_dimensions)
        if source_atlas_dimensions
        else None,
        "selected_metadata_columns": list(selected_metadata_columns),
        "receptor_row_count": len(receptors),
        "unique_cell_count": int(receptors["cell_id"].nunique()),
        "chain_counts": {str(k): int(v) for k, v in receptors["chain"].value_counts().items()},
        "productive_row_count": int(receptors["productive"].fillna(False).sum()),
        "extraction_timestamp": datetime.now(timezone.utc).isoformat(),
        "scrfu_version": __version__,
        "adapter_configuration": adapter_configuration or {},
        "adapter_qc": adapter_qc or {},
        "runtime_versions": {"python": platform.python_version(), "pandas": pd.__version__},
        "files": {
            "receptors": {"name": receptor_path.name, "sha256": file_sha256(receptor_path)},
            "obs_metadata": {"name": metadata_path.name, "sha256": file_sha256(metadata_path)},
        },
    }
    manifest_path = output / "preparation_manifest.json"
    temporary = manifest_path.with_name(manifest_path.name + ".tmp")
    temporary.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(manifest_path)
    return manifest


def validate_receptor_cache(
    cache_dir: PathLike,
    *,
    source_path: PathLike | None = None,
) -> dict[str, Any]:
    cache = Path(cache_dir).expanduser().resolve()
    manifest_path = cache / "preparation_manifest.json"
    errors: list[str] = []
    if not manifest_path.is_file():
        return {"status": "invalid", "structurally_valid": False, "errors": ["manifest missing"]}
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("cache_schema_version") != RECEPTOR_CACHE_SCHEMA_VERSION:
        if (cache / "tcr_ir.tsv.gz").is_file():
            return {
                "status": "legacy_wells_cache",
                "structurally_valid": False,
                "checksum_valid": None,
                "source_status": "not_checked",
                "errors": ["Legacy Wells cache; use migrate_wells_receptor_cache()."],
            }
        errors.append(f"unsupported cache schema version {manifest.get('cache_schema_version')!r}")
    if manifest.get("receptor_schema_version") != RECEPTOR_SCHEMA_VERSION:
        errors.append("unsupported receptor schema version")
    checksum_valid = True
    for key in ("receptors", "obs_metadata"):
        entry = manifest.get("files", {}).get(key, {})
        path = cache / str(entry.get("name", ""))
        if not path.is_file():
            errors.append(f"cache file missing: {path.name}")
            checksum_valid = False
        elif file_sha256(path) != entry.get("sha256"):
            errors.append(f"cache file checksum mismatch: {path.name}")
            checksum_valid = False
    source_status = "unavailable" if source_path is None else "unchanged"
    if source_path is not None:
        source = Path(source_path).expanduser().resolve()
        if not source.is_file():
            source_status = "unavailable"
        else:
            actual = source_fingerprint(source)
            expected = manifest.get("source") or {}
            keys = ("size_bytes", "mtime_ns", "fingerprint_method", "fingerprint")
            if any(actual.get(key) != expected.get(key) for key in keys):
                source_status = "changed"
                errors.append("source fingerprint changed")
    return {
        "status": "valid" if not errors else "invalid",
        "structurally_valid": not any("schema" in error or "missing" in error for error in errors),
        "checksum_valid": checksum_valid,
        "source_status": source_status,
        "errors": errors,
        "manifest": manifest,
    }


def read_receptor_cache(
    cache_dir: PathLike,
    *,
    source_path: PathLike | None = None,
) -> ReceptorCacheData:
    validation = validate_receptor_cache(cache_dir, source_path=source_path)
    if validation["status"] == "legacy_wells_cache":
        raise ValueError(
            "Legacy Wells cache detected. Load it with scrfu.wells.load_wells_receptor_cache() "
            "or migrate it explicitly with scrfu.io.migrate_wells_receptor_cache()."
        )
    if validation["status"] != "valid":
        raise ValueError("Invalid receptor cache: " + "; ".join(validation["errors"]))
    cache = Path(cache_dir).expanduser().resolve()
    manifest = validation["manifest"]
    receptors = pd.read_csv(
        cache / manifest["files"]["receptors"]["name"], sep="\t", dtype={"cell_id": "string"}
    )
    metadata = pd.read_csv(
        cache / manifest["files"]["obs_metadata"]["name"], sep="\t", dtype={"cell_id": "string"}
    )
    validate_receptor_table(receptors, strict=True)
    if len(receptors) != int(manifest["receptor_row_count"]):
        raise ValueError("Receptor cache row count does not match its manifest.")
    metadata = _metadata_for_cache(metadata).set_index("cell_id")
    unknown = sorted(set(receptors["cell_id"].astype(str)).difference(metadata.index.astype(str)))
    if unknown and len(metadata):
        raise ValueError(f"Receptor cache metadata alignment failed for cells: {unknown[:10]}")
    return ReceptorCacheData(receptors, metadata, manifest, validation)


def migrate_wells_receptor_cache(
    legacy_cache_dir: PathLike,
    output_dir: PathLike,
    *,
    force: bool = False,
) -> dict[str, Any]:
    """Explicitly migrate a legacy cache without overwriting it."""
    from .adapters import adapt_wells_tcr_ir
    from .wells import load_wells_receptor_cache

    legacy = load_wells_receptor_cache(legacy_cache_dir)
    result = adapt_wells_tcr_ir(legacy, primary_chain=False)
    return write_receptor_cache(
        output_dir,
        result.receptors,
        result.cell_metadata,
        source_adapter=result.adapter_name,
        source_adapter_version=result.adapter_version,
        source_format="legacy_wells_cache",
        adapter_qc=result.qc,
        adapter_configuration={"primary_chain": False, "migrated_from": str(legacy_cache_dir)},
        source_atlas_dimensions=legacy.atlas_shape,
        force=force,
    )


def export_rfu_matrix(
    adata: object,
    groupby: str,
    path: PathLike,
    normalize: bool = True,
    sep: str = "\t",
) -> pd.DataFrame:
    """
    Export a group-by-RFU matrix to TSV or CSV.
    """
    matrix = aggregate_rfu(adata, groupby=groupby, normalize=normalize)
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(output_path, sep=sep)
    return matrix


__all__ = [
    "H5ADDataFrameInfo",
    "RECEPTOR_CACHE_SCHEMA_VERSION",
    "ReceptorCacheData",
    "UnsupportedH5ADLayout",
    "export_rfu_matrix",
    "file_sha256",
    "migrate_wells_receptor_cache",
    "read_h5ad",
    "read_h5ad_dataframe",
    "read_h5ad_obs",
    "read_h5ad_shape",
    "read_receptor_cache",
    "source_fingerprint",
    "validate_receptor_cache",
    "write_h5ad",
    "write_receptor_cache",
]
