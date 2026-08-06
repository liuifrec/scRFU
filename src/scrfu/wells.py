from __future__ import annotations

import hashlib
import json
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ._version import __version__

WELLS_CACHE_SCHEMA_VERSION = 1
_FINGERPRINT_SAMPLE_BYTES = 1024 * 1024


class UnsupportedWellsH5ADLayout(ValueError):
    """Raised when a Wells H5AD does not use a supported encoded dataframe layout."""


@dataclass
class WellsReceptorData:
    """Expression-free Wells receptor data aligned to selected observations."""

    obs: pd.DataFrame
    tcr_ir: pd.DataFrame
    atlas_shape: tuple[int, int]
    tcr_ir_container: str
    source_h5ad: Path | None = None
    cache_manifest: dict[str, Any] | None = None

    @property
    def obs_names(self) -> pd.Index:
        return self.obs.index

    @property
    def uns(self) -> dict[str, pd.DataFrame]:
        return {"TCR_IR": self.tcr_ir}

    @property
    def obsm(self) -> dict[str, Any]:
        return {}


def _text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _string_list(values: Any) -> list[str]:
    if isinstance(values, (str, bytes)):
        return [_text(values)]
    return [_text(value) for value in values]


def _encoding_type(element: Any) -> str:
    return _text(element.attrs.get("encoding-type", ""))


def _encoding_version(element: Any) -> str:
    return _text(element.attrs.get("encoding-version", ""))


def _validate_dataframe_group(group: Any, *, logical_name: str) -> tuple[str, list[str], int]:
    import h5py

    if not isinstance(group, h5py.Group):
        raise UnsupportedWellsH5ADLayout(
            f"{logical_name} must be an HDF5 group encoded as an AnnData dataframe."
        )
    if _encoding_type(group) != "dataframe":
        seen = _encoding_type(group) or "missing"
        raise UnsupportedWellsH5ADLayout(
            f"{logical_name} has unsupported encoding-type {seen!r}; expected 'dataframe'."
        )
    if _encoding_version(group) != "0.2.0":
        seen = _encoding_version(group) or "missing"
        raise UnsupportedWellsH5ADLayout(
            f"{logical_name} has unsupported dataframe encoding-version {seen!r}; expected '0.2.0'."
        )
    index_key = _text(group.attrs.get("_index", ""))
    if not index_key or index_key not in group:
        raise UnsupportedWellsH5ADLayout(
            f"{logical_name} dataframe does not declare a readable '_index' dataset."
        )
    index = group[index_key]
    if not isinstance(index, h5py.Dataset) or index.ndim != 1:
        raise UnsupportedWellsH5ADLayout(
            f"{logical_name} dataframe index must be a one-dimensional dataset."
        )
    if _encoding_type(index) not in {"array", "string-array"}:
        seen = _encoding_type(index) or "missing"
        raise UnsupportedWellsH5ADLayout(
            f"{logical_name} dataframe index has unsupported encoding-type {seen!r}."
        )
    if "column-order" not in group.attrs:
        raise UnsupportedWellsH5ADLayout(
            f"{logical_name} dataframe is missing its 'column-order' attribute."
        )
    columns = _string_list(group.attrs["column-order"])
    missing = [column for column in columns if column not in group]
    if missing:
        raise UnsupportedWellsH5ADLayout(
            f"{logical_name} dataframe declares missing columns: {missing}."
        )
    return index_key, columns, int(index.shape[0])


def _take(dataset: Any, positions: np.ndarray) -> np.ndarray:
    """Read arbitrary unique rows while satisfying h5py's sorted-index requirement."""
    if len(positions) == 0:
        return np.asarray(dataset[:0])
    order = np.argsort(positions, kind="stable")
    sorted_positions = positions[order]
    if len(np.unique(sorted_positions)) != len(sorted_positions):
        raise ValueError("Selected observation positions must be unique.")
    values = np.asarray(dataset[sorted_positions])
    inverse = np.empty_like(order)
    inverse[order] = np.arange(len(order))
    return values[inverse]


def _decode_array(values: np.ndarray) -> np.ndarray:
    if values.dtype.kind in {"S", "O", "U"}:
        return np.asarray(
            [value.decode("utf-8") if isinstance(value, bytes) else value for value in values],
            dtype=object,
        )
    return values


def _read_encoded_column(element: Any, positions: np.ndarray, *, logical_name: str) -> Any:
    import h5py

    if isinstance(element, h5py.Dataset):
        encoding = _encoding_type(element)
        if encoding not in {"array", "string-array"}:
            seen = encoding or "missing"
            raise UnsupportedWellsH5ADLayout(
                f"{logical_name} has unsupported dataset encoding-type {seen!r}."
            )
        if element.ndim != 1:
            raise UnsupportedWellsH5ADLayout(
                f"{logical_name} must be one-dimensional; found shape {element.shape}."
            )
        return _decode_array(_take(element, positions))

    if not isinstance(element, h5py.Group):
        raise UnsupportedWellsH5ADLayout(f"{logical_name} is not a supported HDF5 element.")

    encoding = _encoding_type(element)
    if encoding == "categorical":
        if _encoding_version(element) != "0.2.0":
            seen = _encoding_version(element) or "missing"
            raise UnsupportedWellsH5ADLayout(
                f"{logical_name} has unsupported categorical encoding-version {seen!r}."
            )
        if "codes" not in element or "categories" not in element:
            raise UnsupportedWellsH5ADLayout(
                f"{logical_name} categorical encoding requires 'codes' and 'categories'."
            )
        codes = _take(element["codes"], positions).astype(np.int64, copy=False)
        categories = _decode_array(np.asarray(element["categories"][:]))
        if np.any(codes < -1) or np.any(codes >= len(categories)):
            raise UnsupportedWellsH5ADLayout(
                f"{logical_name} contains categorical codes outside the category range."
            )
        ordered = bool(element.attrs.get("ordered", False))
        return pd.Categorical.from_codes(codes, categories=categories, ordered=ordered)

    if encoding in {"nullable-integer", "nullable-boolean", "nullable-string-array"}:
        if "values" not in element or "mask" not in element:
            raise UnsupportedWellsH5ADLayout(
                f"{logical_name} {encoding!r} encoding requires 'values' and 'mask'."
            )
        values = _decode_array(_take(element["values"], positions))
        mask = _take(element["mask"], positions).astype(bool, copy=False)
        series = pd.Series(values)
        series.loc[mask] = pd.NA
        if encoding == "nullable-integer":
            return pd.array(series, dtype="Int64")
        if encoding == "nullable-boolean":
            return pd.array(series, dtype="boolean")
        return pd.array(series, dtype="string")

    seen = encoding or "missing"
    raise UnsupportedWellsH5ADLayout(f"{logical_name} has unsupported encoding-type {seen!r}.")


def _read_dataframe(
    group: Any,
    *,
    logical_name: str,
    positions: np.ndarray,
    columns: Sequence[str] | None,
) -> pd.DataFrame:
    index_key, available, row_count = _validate_dataframe_group(group, logical_name=logical_name)
    if np.any(positions < 0) or np.any(positions >= row_count):
        raise IndexError(f"Selected positions fall outside {logical_name} with {row_count} rows.")
    selected_columns = available if columns is None else list(columns)
    missing = [column for column in selected_columns if column not in available]
    if missing:
        raise ValueError(
            f"Requested columns are missing from {logical_name}: {missing}. "
            f"Available columns: {available}"
        )
    index_values = _decode_array(_take(group[index_key], positions))
    index = pd.Index([_text(value) for value in index_values], name=None)
    data = {
        column: _read_encoded_column(
            group[column], positions, logical_name=f"{logical_name}[{column!r}]"
        )
        for column in selected_columns
    }
    return pd.DataFrame(data, index=index)


def _atlas_shape(handle: Any, obs_count: int) -> tuple[int, int]:
    import h5py

    if "X" not in handle:
        if "var" not in handle:
            raise UnsupportedWellsH5ADLayout(
                "H5AD root has neither X nor the var dataframe; atlas dimensions are unknown."
            )
        _, _, var_count = _validate_dataframe_group(handle["var"], logical_name="var")
        return obs_count, var_count
    x = handle["X"]
    if isinstance(x, h5py.Dataset):
        shape = tuple(map(int, x.shape))
    elif isinstance(x, h5py.Group) and "shape" in x.attrs:
        shape = tuple(map(int, x.attrs["shape"]))
    else:
        raise UnsupportedWellsH5ADLayout(
            "H5AD 'X' has no supported dataset shape or sparse-group 'shape' attribute."
        )
    if len(shape) != 2 or shape[0] != obs_count:
        raise UnsupportedWellsH5ADLayout(
            f"H5AD X shape {shape} is inconsistent with {obs_count} observation rows."
        )
    return shape


def _select_positions(
    obs_group: Any,
    *,
    obs_count: int,
    max_cells: int | None,
    selected_obs_names: Sequence[str] | None,
) -> np.ndarray:
    if max_cells is not None and selected_obs_names is not None:
        raise ValueError("Pass either max_cells or selected_obs_names, not both.")
    if max_cells is not None:
        if isinstance(max_cells, bool) or not isinstance(max_cells, (int, np.integer)):
            raise ValueError("max_cells must be a positive integer.")
        if max_cells <= 0:
            raise ValueError("max_cells must be a positive integer.")
        return np.arange(min(int(max_cells), obs_count), dtype=np.int64)
    if selected_obs_names is None:
        return np.arange(obs_count, dtype=np.int64)

    index_key, _, _ = _validate_dataframe_group(obs_group, logical_name="obs")
    all_names = [_text(value) for value in _decode_array(np.asarray(obs_group[index_key][:]))]
    lookup = {name: position for position, name in enumerate(all_names)}
    if len(lookup) != len(all_names):
        raise UnsupportedWellsH5ADLayout("obs index contains duplicate observation names.")
    requested = list(map(str, selected_obs_names))
    missing = [name for name in requested if name not in lookup]
    if missing:
        raise KeyError(f"Selected observation names are absent from the H5AD: {missing[:10]}")
    return np.asarray([lookup[name] for name in requested], dtype=np.int64)


def read_wells_receptors_h5ad(
    source_h5ad: str | Path,
    *,
    obs_columns: Sequence[str] = (),
    max_cells: int | None = None,
    selected_obs_names: Sequence[str] | None = None,
) -> WellsReceptorData:
    """Read Wells receptors and selected metadata without reading expression matrices.

    The reader validates AnnData's encoded-dataframe schema at runtime. It reads
    ``TCR_IR`` from ``uns`` (preferred) or ``obsm`` only after confirming that the
    selected element is a row-aligned encoded dataframe.
    """
    import h5py

    path = Path(source_h5ad).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Wells H5AD not found: {path}")
    with h5py.File(path, "r") as handle:
        if _encoding_type(handle) != "anndata":
            seen = _encoding_type(handle) or "missing"
            raise UnsupportedWellsH5ADLayout(
                f"H5AD root has encoding-type {seen!r}; expected 'anndata'."
            )
        if "obs" not in handle:
            raise UnsupportedWellsH5ADLayout("H5AD root is missing the obs dataframe.")
        obs_group = handle["obs"]
        _, _, obs_count = _validate_dataframe_group(obs_group, logical_name="obs")
        atlas_shape = _atlas_shape(handle, obs_count)

        candidates = [
            ("uns", handle.get("uns/TCR_IR")),
            ("obsm", handle.get("obsm/TCR_IR")),
        ]
        container, tcr_group = next(
            ((name, group) for name, group in candidates if group is not None), ("", None)
        )
        if tcr_group is None:
            raise UnsupportedWellsH5ADLayout(
                "H5AD contains neither uns['TCR_IR'] nor obsm['TCR_IR']."
            )
        _, _, tcr_count = _validate_dataframe_group(
            tcr_group, logical_name=f"{container}['TCR_IR']"
        )
        if tcr_count != obs_count:
            raise UnsupportedWellsH5ADLayout(
                f"{container}['TCR_IR'] has {tcr_count} rows but obs has {obs_count}; "
                "row alignment is unsupported."
            )

    # Interpretation remains Wells-specific, while byte-level dataframe reading is generic.
    from .io import UnsupportedH5ADLayout, read_h5ad_dataframe, read_h5ad_obs

    try:
        obs = read_h5ad_obs(
            path,
            columns=obs_columns,
            selected_names=selected_obs_names,
            max_rows=max_cells,
        )
        tcr_ir = read_h5ad_dataframe(
            path,
            location=container,
            key="TCR_IR",
            selected_names=obs.index.astype(str).tolist(),
        )
    except UnsupportedH5ADLayout as exc:
        raise UnsupportedWellsH5ADLayout(str(exc)) from exc
    if not tcr_ir.index.equals(obs.index):
        raise UnsupportedWellsH5ADLayout(
            f"{container}['TCR_IR'] index does not exactly match the selected obs index."
        )
    return WellsReceptorData(
        obs=obs,
        tcr_ir=tcr_ir,
        atlas_shape=atlas_shape,
        tcr_ir_container=container,
        source_h5ad=path,
    )


def source_fingerprint(source_h5ad: str | Path) -> dict[str, Any]:
    """Return a fast content fingerprint using fixed first/middle/last samples."""
    path = Path(source_h5ad).expanduser().resolve()
    stat = path.stat()
    size = stat.st_size
    offsets = sorted(
        {
            0,
            max(0, size // 2 - _FINGERPRINT_SAMPLE_BYTES // 2),
            max(0, size - _FINGERPRINT_SAMPLE_BYTES),
        }
    )
    digest = hashlib.sha256()
    digest.update(f"size:{size}\n".encode())
    with path.open("rb") as stream:
        for offset in offsets:
            stream.seek(offset)
            block = stream.read(_FINGERPRINT_SAMPLE_BYTES)
            digest.update(f"offset:{offset}:length:{len(block)}\n".encode())
            digest.update(block)
    return {
        "path": str(path),
        "size_bytes": size,
        "mtime_ns": stat.st_mtime_ns,
        "fingerprint_algorithm": "sha256-size-first-middle-last-1MiB-v1",
        "fingerprint": digest.hexdigest(),
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def prepare_wells_receptor_cache(
    source_h5ad: str | Path,
    output_dir: str | Path,
    *,
    obs_columns: Sequence[str] = (),
    max_cells: int | None = None,
) -> dict[str, Any]:
    """Create a compact, expression-free Wells receptor cache."""
    source = Path(source_h5ad).expanduser().resolve()
    output = Path(output_dir).expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    data = read_wells_receptors_h5ad(
        source,
        obs_columns=obs_columns,
        max_cells=max_cells,
    )

    tcr_path = output / "tcr_ir.tsv.gz"
    obs_path = output / "obs_metadata.tsv.gz"
    data.tcr_ir.rename_axis("cell_id").reset_index().to_csv(
        tcr_path, sep="\t", index=False, compression="gzip"
    )
    data.obs.rename_axis("cell_id").reset_index().to_csv(
        obs_path, sep="\t", index=False, compression="gzip"
    )
    manifest = {
        "cache_schema_version": WELLS_CACHE_SCHEMA_VERSION,
        "source_h5ad": source_fingerprint(source),
        "source_atlas_dimensions": list(data.atlas_shape),
        "selected_cell_count": len(data.obs),
        "selected_metadata_columns": list(obs_columns),
        "tcr_ir_container": data.tcr_ir_container,
        "tcr_row_count": len(data.tcr_ir),
        "extraction_timestamp": datetime.now(timezone.utc).isoformat(),
        "scrfu_version": __version__,
        "files": {
            "tcr_ir": {"name": tcr_path.name, "sha256": _sha256(tcr_path)},
            "obs_metadata": {"name": obs_path.name, "sha256": _sha256(obs_path)},
        },
    }
    _atomic_json(output / "preparation_manifest.json", manifest)
    return manifest


def load_wells_receptor_cache(
    cache_dir: str | Path,
    *,
    source_h5ad: str | Path | None = None,
) -> WellsReceptorData:
    """Load and validate a compact Wells receptor cache."""
    cache = Path(cache_dir).expanduser().resolve()
    manifest_path = cache / "preparation_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Wells cache manifest not found: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("cache_schema_version") != WELLS_CACHE_SCHEMA_VERSION:
        raise ValueError(
            "Unsupported Wells cache schema version: "
            f"{manifest.get('cache_schema_version')!r}; expected {WELLS_CACHE_SCHEMA_VERSION}."
        )
    if source_h5ad is not None:
        actual = source_fingerprint(source_h5ad)
        expected = manifest.get("source_h5ad", {})
        keys = ("size_bytes", "mtime_ns", "fingerprint_algorithm", "fingerprint")
        mismatches = [key for key in keys if actual.get(key) != expected.get(key)]
        if mismatches:
            raise ValueError(
                "Wells receptor cache is stale for the supplied source H5AD; "
                f"fingerprint fields changed: {mismatches}. Re-run prepare-wells."
            )

    loaded: dict[str, pd.DataFrame] = {}
    for key in ("tcr_ir", "obs_metadata"):
        entry = manifest.get("files", {}).get(key, {})
        path = cache / str(entry.get("name", ""))
        if not path.is_file():
            raise ValueError(f"Wells cache file is missing: {path}")
        if _sha256(path) != entry.get("sha256"):
            raise ValueError(f"Wells cache file checksum mismatch: {path.name}")
        frame = pd.read_csv(path, sep="\t", dtype={"cell_id": "string"})
        if "cell_id" not in frame:
            raise ValueError(f"Wells cache file {path.name} is missing 'cell_id'.")
        if frame["cell_id"].isna().any() or frame["cell_id"].duplicated().any():
            raise ValueError(f"Wells cache file {path.name} has invalid cell identifiers.")
        loaded[key] = frame.set_index("cell_id")

    obs = loaded["obs_metadata"]
    tcr_ir = loaded["tcr_ir"]
    if not tcr_ir.index.equals(obs.index):
        raise ValueError("Wells cache TCR_IR rows do not exactly align with obs metadata rows.")
    expected_rows = int(manifest.get("tcr_row_count", -1))
    if len(tcr_ir) != expected_rows or len(obs) != int(manifest.get("selected_cell_count", -1)):
        raise ValueError("Wells cache row counts do not match the preparation manifest.")
    shape = manifest.get("source_atlas_dimensions", [])
    if not isinstance(shape, list) or len(shape) != 2:
        raise ValueError("Wells cache manifest has invalid source_atlas_dimensions.")
    return WellsReceptorData(
        obs=obs,
        tcr_ir=tcr_ir,
        atlas_shape=(int(shape[0]), int(shape[1])),
        tcr_ir_container=str(manifest.get("tcr_ir_container", "uns")),
        source_h5ad=Path(source_h5ad).resolve() if source_h5ad is not None else None,
        cache_manifest=manifest,
    )


__all__ = [
    "UnsupportedWellsH5ADLayout",
    "WELLS_CACHE_SCHEMA_VERSION",
    "WellsReceptorData",
    "load_wells_receptor_cache",
    "prepare_wells_receptor_cache",
    "read_wells_receptors_h5ad",
    "source_fingerprint",
]
