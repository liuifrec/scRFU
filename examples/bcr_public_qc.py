"""Run experimental BCR canonicalization and feature QC on public Cell Ranger tables."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from scrfu import __version__, bcr


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _library(path: Path) -> str:
    return re.sub(r"_all_contig_annotations\.csv\.gz$", "", path.name)


def _sample_key(library: str) -> str:
    return re.sub(r"^GSM\d+_", "", library)


def _read_inputs(directory: Path, pattern: str) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    tables: list[pd.DataFrame] = []
    artifacts: list[dict[str, Any]] = []
    for path in sorted(directory.glob(pattern)):
        frame = pd.read_csv(path, low_memory=False)
        library = _library(path)
        frame["library_id"] = library
        for column in ("barcode", "contig_id", "raw_clonotype_id", "raw_consensus_id"):
            if column in frame:
                values = frame[column].astype("string")
                frame[column] = values.where(values.isna(), library + "::" + values)
        tables.append(frame)
        artifacts.append(
            {"filename": path.name, "size_bytes": path.stat().st_size, "sha256": _sha256(path)}
        )
    if not tables:
        raise FileNotFoundError(f"No public BCR files match {directory / pattern}")
    return pd.concat(tables, ignore_index=True), artifacts


def _demux_qc(directory: Path | None, cells: pd.Series) -> dict[str, Any]:
    if directory is None:
        return {"configured": False}
    cell_keys = {re.sub(r"^GSM\d+_", "", value) for value in cells.astype(str)}
    assigned: set[str] = set()
    donor_labels: set[str] = set()
    files = sorted(directory.glob("*.best"))
    for path in files:
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        sample = path.name.removesuffix("_demux.best")
        keys = sample + "::" + frame["BARCODE"].astype(str)
        single = frame["DROPLET.TYPE"].astype(str).eq("SNG") & keys.isin(cell_keys)
        assigned.update(keys.loc[single].astype(str))
        donor_labels.update(frame.loc[single, "SNG.BEST.GUESS"].dropna().astype(str))
    return {
        "configured": True,
        "file_count": len(files),
        "single_donor_cell_count": len(assigned),
        "single_donor_cell_fraction": len(assigned) / max(1, len(cell_keys)),
        "donor_count": len(donor_labels),
        "donor_labels_exported": False,
    }


def run_qc(
    *,
    dataset_label: str,
    input_dir: Path,
    pattern: str,
    output_dir: Path,
    demux_dir: Path | None,
) -> dict[str, Any]:
    raw, inputs = _read_inputs(input_dir, pattern)
    input_rows = len(raw)
    is_cell = raw.get("is_cell", pd.Series(True, index=raw.index)).eq(True)
    confidence = raw.get("high_confidence", pd.Series(True, index=raw.index)).eq(True)
    productive = raw.get("productive", pd.Series(True, index=raw.index)).eq(True)
    selected = raw.loc[is_cell & confidence & productive].copy()
    prepared = bcr.prepare_bcr_table(selected, source_label=dataset_label)
    feature_result = bcr.bcr_feature_matrix(prepared.receptors, pairs=prepared.pairs)
    receptors = prepared.receptors
    pairs = prepared.pairs
    cdr3 = receptors["cdr3aa"].astype("string")
    duplicate_rows = int(cdr3.notna().sum() - cdr3.dropna().nunique())
    qc = {
        "dataset_label": dataset_label,
        "input_rows": input_rows,
        "is_cell_high_confidence_rows": int((is_cell & confidence).sum()),
        "productive_rows": len(receptors),
        "chain_counts": {
            str(key): int(value) for key, value in receptors["chain"].value_counts().items()
        },
        "cell_count": int(receptors["cell_id"].nunique()),
        "pair_status_counts": {
            str(key): int(value) for key, value in pairs["pair_status"].value_counts().items()
        },
        "v_complete_fraction": float(receptors["v_call"].notna().mean()),
        "j_complete_fraction": float(receptors["j_call"].notna().mean()),
        "constant_complete_fraction": float(receptors["c_call"].notna().mean()),
        "isotype_complete_fraction": float(receptors["isotype"].notna().mean()),
        "mutation_frequency_complete_fraction": float(
            receptors["mutation_frequency"].notna().mean()
        ),
        "germline_identity_complete_fraction": float(receptors["germline_identity"].notna().mean()),
        "clonotype_complete_fraction": float(receptors["clonotype_id"].notna().mean()),
        "clonal_family_complete_fraction": float(receptors["clonal_family_id"].notna().mean()),
        "duplicate_exact_cdr3_rows": duplicate_rows,
        "duplicate_exact_cdr3_fraction": duplicate_rows / max(1, int(cdr3.notna().sum())),
        "malformed_or_missing_cdr3_rows": int((~cdr3.str.match(r"^C[A-Z*]+$", na=False)).sum()),
        "cells_with_multiple_heavy_candidates": int(pairs["heavy_candidate_count"].gt(1).sum()),
        "cells_with_multiple_light_candidates": int(pairs["light_candidate_count"].gt(1).sum()),
        "source_sequence_ids_traceable": bool(
            receptors["sequence_id"].astype(str).str.contains("::", regex=False).all()
        ),
        "one_feature_row_per_cell": bool(
            len(feature_result.features) == receptors["cell_id"].nunique()
            and feature_result.features["cell_id"].is_unique
        ),
        "demultiplexing": _demux_qc(demux_dir, receptors["cell_id"]),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    feature_path = output_dir / "bcr_feature_matrix.tsv.gz"
    missingness_path = output_dir / "bcr_feature_missingness.tsv"
    qc_path = output_dir / "bcr_qc.json"
    feature_result.features.to_csv(feature_path, sep="\t", index=False)
    feature_result.missingness.to_csv(missingness_path, sep="\t", index=False)
    qc_path.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n")
    outputs = [feature_path, missingness_path, qc_path]
    manifest = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset_label": dataset_label,
        "scrfu_version": __version__,
        "experimental": True,
        "input_artifacts": inputs,
        "parameters": {
            "pattern": pattern,
            "is_cell": True,
            "high_confidence": True,
            "productive": True,
            "outcome_labels_used": False,
            "functional_reference_built": False,
        },
        "qc": qc,
        "outputs": {
            path.name: {
                "size_bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
            for path in outputs
        },
    }
    (output_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-label", required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--pattern", default="*contig_annotations.csv.gz")
    parser.add_argument("--demux-dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    report = run_qc(
        dataset_label=args.dataset_label,
        input_dir=args.input_dir,
        pattern=args.pattern,
        output_dir=args.output_dir,
        demux_dir=args.demux_dir,
    )
    print(json.dumps(report["qc"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
