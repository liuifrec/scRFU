from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import pandas as pd

from .repertoire import analysis_frame, repertoire_metrics


@dataclass(frozen=True)
class ComparatorRepresentation:
    name: str
    matrix: pd.DataFrame
    sample_metadata: pd.DataFrame
    parameters: dict[str, Any]


Comparator = Callable[..., ComparatorRepresentation]
_EXTERNAL_COMPARATORS: dict[str, Comparator] = {}


def register_comparator(name: str, comparator: Comparator, *, replace: bool = False) -> None:
    """Register an explicit optional comparator without adding a package dependency."""
    normalized = str(name).strip().lower()
    if not normalized or not callable(comparator):
        raise ValueError("Comparator name and callable are required.")
    if normalized in _EXTERNAL_COMPARATORS and not replace:
        raise ValueError(f"Comparator {normalized!r} is already registered.")
    _EXTERNAL_COMPARATORS[normalized] = comparator


def list_comparators() -> list[str]:
    builtins = [
        "cdr3_length",
        "clonotype",
        "dominant_clone_fraction",
        "edit_distance",
        "exact_cdr3",
        "j_gene",
        "shannon",
        "simpson",
        "v_gene",
    ]
    return sorted(set(builtins) | set(_EXTERNAL_COMPARATORS))


def _sample_metadata(frame: pd.DataFrame, sample_key: str) -> pd.DataFrame:
    return pd.DataFrame(
        index=pd.Index(sorted(frame[sample_key].dropna().unique(), key=str), name=sample_key)
    )


def _feature_matrix(
    frame: pd.DataFrame, sample_key: str, feature: pd.Series, cell_col: str
) -> pd.DataFrame:
    source = pd.DataFrame(
        {sample_key: frame[sample_key], cell_col: frame[cell_col], "feature": feature},
        index=frame.index,
    )
    source = source.dropna().drop_duplicates([sample_key, cell_col, "feature"])
    matrix = source.groupby([sample_key, "feature"], observed=True).size().unstack(fill_value=0)
    matrix = matrix.reindex(sorted(matrix.index, key=str))
    matrix = matrix.reindex(sorted(matrix.columns, key=str), axis=1)
    matrix.columns = matrix.columns.astype(str)
    return matrix


def _edit_distance(left: str, right: str) -> int:
    previous = list(range(len(right) + 1))
    for left_index, left_character in enumerate(left, start=1):
        current = [left_index]
        for right_index, right_character in enumerate(right, start=1):
            current.append(
                min(
                    current[-1] + 1,
                    previous[right_index] + 1,
                    previous[right_index - 1] + (left_character != right_character),
                )
            )
        previous = current
    return previous[-1]


def _edit_distance_labels(sequences: pd.Series, threshold: int) -> pd.Series:
    unique = sorted(sequences.dropna().astype(str).unique())
    parent = {sequence: sequence for sequence in unique}

    def find(sequence: str) -> str:
        while parent[sequence] != sequence:
            parent[sequence] = parent[parent[sequence]]
            sequence = parent[sequence]
        return sequence

    for left_index, left in enumerate(unique):
        for right in unique[left_index + 1 :]:
            if (
                abs(len(left) - len(right)) <= threshold
                and _edit_distance(left, right) <= threshold
            ):
                left_root, right_root = find(left), find(right)
                if left_root != right_root:
                    parent[max(left_root, right_root)] = min(left_root, right_root)
    roots = {sequence: find(sequence) for sequence in unique}
    root_ids = {root: index for index, root in enumerate(sorted(set(roots.values())))}
    labels = {sequence: f"edit_cluster_{root_ids[root]:06d}" for sequence, root in roots.items()}
    return sequences.astype("string").map(labels)


def repertoire_representation(
    data: Any,
    *,
    sample_key: str,
    method: str = "exact_cdr3",
    cdr3_col: str = "cdr3aa",
    clonotype_col: str = "clonotype_id",
    v_col: str = "v_call",
    j_col: str = "j_call",
    cell_col: str = "cell_id",
    max_edit_distance: int = 1,
) -> ComparatorRepresentation:
    """Build a sample-level conventional comparator for shared evaluation interfaces."""
    normalized = method.strip().lower()
    if normalized in _EXTERNAL_COMPARATORS:
        return _EXTERNAL_COMPARATORS[normalized](data=data, sample_key=sample_key)
    frame = analysis_frame(data)
    if sample_key not in frame:
        raise ValueError(f"Comparator input is missing sample key {sample_key!r}.")
    if cell_col not in frame:
        raise ValueError(f"Comparator input is missing cell key {cell_col!r}.")
    metadata = _sample_metadata(frame, sample_key)
    feature_columns = {
        "exact_cdr3": cdr3_col,
        "clonotype": clonotype_col,
        "v_gene": v_col,
        "j_gene": j_col,
    }
    if normalized in feature_columns:
        column = feature_columns[normalized]
        if column not in frame:
            raise ValueError(f"Comparator {normalized!r} requires {column!r}.")
        matrix = _feature_matrix(frame, sample_key, frame[column], cell_col)
    elif normalized == "cdr3_length":
        if cdr3_col not in frame:
            raise ValueError(f"Comparator requires {cdr3_col!r}.")
        matrix = _feature_matrix(
            frame, sample_key, frame[cdr3_col].astype("string").str.len(), cell_col
        )
    elif normalized == "edit_distance":
        if isinstance(max_edit_distance, bool) or not isinstance(max_edit_distance, int):
            raise TypeError("max_edit_distance must be a non-negative integer.")
        if max_edit_distance < 0 or cdr3_col not in frame:
            raise ValueError("Edit-distance comparison requires CDR3 and a non-negative threshold.")
        matrix = _feature_matrix(
            frame,
            sample_key,
            _edit_distance_labels(frame[cdr3_col], max_edit_distance),
            cell_col,
        )
    elif normalized in {"shannon", "simpson", "dominant_clone_fraction"}:
        metric_name = {
            "shannon": "shannon_entropy",
            "simpson": "simpson_diversity",
            "dominant_clone_fraction": "dominant_clonotype_fraction",
        }[normalized]
        summary = repertoire_metrics(
            frame,
            groupby=sample_key,
            weighting="cell",
            cell_col=cell_col,
            cdr3_col=cdr3_col,
            clonotype_col=clonotype_col,
        )
        if metric_name not in summary:
            if normalized == "dominant_clone_fraction":
                raise ValueError("Dominant-clone comparison requires true clonotype identifiers.")
            raise ValueError(f"Repertoire metrics lack {metric_name!r}.")
        matrix = summary.set_index(sample_key)[[metric_name]].rename(
            columns={metric_name: normalized}
        )
    else:
        raise ValueError(f"Unknown comparator {method!r}. Available: {list_comparators()}")
    matrix = matrix.reindex(metadata.index, fill_value=0).astype(float)
    return ComparatorRepresentation(
        normalized,
        matrix,
        metadata,
        {
            "sample_key": sample_key,
            "method": normalized,
            "unit": "biological_sample",
            "external": False,
            "max_edit_distance": max_edit_distance if normalized == "edit_distance" else None,
        },
    )


__all__ = [
    "ComparatorRepresentation",
    "list_comparators",
    "register_comparator",
    "repertoire_representation",
]
