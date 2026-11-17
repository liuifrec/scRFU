"""Prepare local RfuWAS Supplementary Data 1 as GRCh38 RFU-QTL evidence.

XLSX inputs select sheet 'Data 1', with headers on row 2 (optional openpyxl
reader). TSV/CSV inputs must be exports of Data 1 with the column header first.
No download, effect-allele inference, or RFU label offset is performed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from decimal import Decimal, InvalidOperation
from pathlib import Path

import pandas as pd

from scrfu.tl import normalize_regulatory_variants

SOURCE = "RfuWAS Supplementary Data 1"
URL = "https://doi.org/10.1038/s42003-024-07010-x"
REQUIRED = ("SNP", "RFU", "beta", "t.stat", "p.value")


def _label(value: object) -> str:
    try:
        number = Decimal(str(value))
    except InvalidOperation as exc:
        raise ValueError("RFU labels must be positive one-based integers") from exc
    if not number.is_finite() or number <= 0 or number != number.to_integral_value():
        raise ValueError("RFU labels must be positive one-based integers")
    return str(int(number))


def prepare_rfuwas_data1(table: pd.DataFrame, *, release: str) -> pd.DataFrame:
    """Validate Data 1, preserving original columns and published RFU identities.

    REF/ALT describe variant identity only. Effect allele remains missing and
    beta is preserved without flipping. SE is not reconstructed from t.stat.
    """
    if not release.strip():
        raise ValueError("An explicit nonempty release is required")
    if not table.columns.is_unique:
        raise ValueError("Duplicate input column names")
    missing = set(REQUIRED) - set(table.columns)
    if missing:
        raise ValueError(f"Expected RfuWAS Data 1 columns; missing {sorted(missing)}")
    generated = {
        "chromosome",
        "position",
        "ref",
        "alt",
        "variant_id",
        "genome_build",
        "rfu_label",
        "pvalue",
        "source",
        "release",
        "url",
        "effect_allele",
        "effect_allele_status",
    }
    if generated.intersection(table.columns):
        raise ValueError("Use the original Data 1 schema, not an already adapted table")
    frame = table.copy()
    parts = (
        frame["SNP"]
        .astype("string")
        .str.extract(
            r"^(?P<chromosome>(?:chr)?(?:[0-9]+|X|Y|MT|M))_(?P<position>[0-9]+)_(?P<ref>[ACGTacgt]+)_(?P<alt>[ACGTacgt]+)$"
        )
    )
    if parts.isna().any(axis=None):
        raise ValueError("SNP must have the form CHR_POS_REF_ALT with biallelic A/C/G/T sequences")
    for column in parts:
        frame[column] = parts[column]
    frame["variant_id"] = parts.astype(str).agg(":".join, axis=1)
    frame["rfu_label"] = frame["RFU"].map(_label)
    frame["pvalue"] = frame["p.value"]
    frame["genome_build"] = "GRCh38"
    frame["source"] = SOURCE
    frame["release"] = release.strip()
    frame["url"] = URL
    frame["effect_allele"] = pd.NA
    frame["effect_allele_status"] = "unresolved_in_published_data1"
    normalized = normalize_regulatory_variants(frame, layer="rfu_qtl")
    normalized["variant_id"] = normalized["variant_key"]
    return normalized


def read_data1(path: Path, *, format: str, nrows: int | None = None) -> pd.DataFrame:
    """Read the verified worksheet/header, or an explicitly selected Data 1 export."""
    if nrows is not None and nrows < 1:
        raise ValueError("nrows must be positive")
    if format == "xlsx":
        try:
            return pd.read_excel(
                path, sheet_name="Data 1", header=1, dtype=object, nrows=nrows, engine="openpyxl"
            )
        except ImportError as exc:
            raise ImportError(
                "XLSX reading requires openpyxl in your analysis environment; alternatively export Data 1 as TSV/CSV."
            ) from exc
    if format not in {"tsv", "csv"}:
        raise ValueError("format must be xlsx, tsv, or csv")
    return pd.read_csv(path, sep="\t" if format == "tsv" else ",", dtype=object, nrows=nrows)


def run(
    input_path: Path, outdir: Path, *, format: str, release: str, nrows: int | None = None
) -> None:
    """Write canonical TSV and deterministic provenance without overwriting files."""
    outputs = [outdir / "rfu_qtl.tsv", outdir / "provenance.json"]
    if any(path.exists() for path in outputs):
        raise FileExistsError("Output files already exist; choose a new output directory")
    raw = read_data1(input_path, format=format, nrows=nrows)
    canonical = prepare_rfuwas_data1(raw, release=release)
    digest = hashlib.sha256()
    with input_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    canonical["sha256"] = digest.hexdigest()
    provenance = {
        "schema_version": 1,
        "source": SOURCE,
        "release": release.strip(),
        "url": URL,
        "input_path": str(input_path.resolve()),
        "input_sha256": digest.hexdigest(),
        "input_format": format,
        "input_object": "Data 1",
        "genome_build": "GRCh38",
        "worksheet": "Data 1" if format == "xlsx" else None,
        "header_row": 2 if format == "xlsx" else 1,
        "nrows_limit": nrows,
        "rows_read": len(raw),
        "rows_written": len(canonical),
        "effect_allele_status": "unresolved_in_published_data1",
        "beta_sign_changed": False,
        "rfu_label_offset": 0,
        "significance_filter": None,
        "input_columns": list(raw.columns),
    }
    outdir.mkdir(parents=True, exist_ok=True)
    canonical.to_csv(outputs[0], sep="\t", index=False)
    outputs[1].write_text(json.dumps(provenance, indent=2, sort_keys=True, allow_nan=False) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--format", required=True, choices=("xlsx", "tsv", "csv"))
    parser.add_argument("--release", required=True)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument(
        "--nrows", type=int, help="Optional positive row limit for a small smoke test"
    )
    args = parser.parse_args()
    run(args.input, args.outdir, format=args.format, release=args.release, nrows=args.nrows)


if __name__ == "__main__":
    main()
