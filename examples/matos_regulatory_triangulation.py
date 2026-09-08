"""Triangulate local RFU-QTL and Matos CD4 QTL exports using a JSON manifest.

See docs/regulatory_triangulation.md for formats, column maps, and a manifest.
Only pandas-readable TSV/CSV (optionally gzip) inputs are read; no downloads.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

from scrfu.tl import normalize_regulatory_variants, regulatory_triangulation

MATOS_URL = "https://doi.org/10.64898/2026.01.27.26344979"
_TARGETS = {"rfu_qtl": "rfu_label", "eqtl": "gene", "caqtl": "peak", "gwas": "trait"}


def adapt_matos_table(
    table: pd.DataFrame,
    *,
    layer: str,
    source: str,
    release: str,
    genome_build: str,
    format: str,
    context: str | None = None,
    allele_order: str | None = None,
    variant_column: str | None = None,
    effect_allele_column: str | None = None,
    column_map: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Adapt an explicit canonical, tensorqtl, or susie export without guessing alleles.

    Custom CHR:POS[b38]A,G IDs require allele_order='ref-alt' or 'alt-ref'.
    The study's ID-generation command and parsing comments differ; verify against
    the release VCF. Effect direction stays unknown without effect_allele_column
    or a canonical effect_allele. Fine-mapping prefers original variant_id.ss.
    """
    if layer not in _TARGETS or format not in {"canonical", "tensorqtl", "susie"}:
        raise ValueError("Unsupported layer or format")
    if format != "canonical" and layer not in {"eqtl", "caqtl"}:
        raise ValueError("Matos formats are only supported for eqtl/caqtl")
    frame = table.rename(columns=column_map or {}).copy()
    if not frame.columns.is_unique:
        raise ValueError("column_map creates duplicate columns")
    if format != "canonical":
        aliases = {
            "phenotype_id": _TARGETS[layer],
            "slope": "beta",
            "slope_se": "se",
            "pval_nominal": "pvalue",
            "af": "allele_frequency",
            "PIP": "pip",
            "cs": "credible_set_id",
            "region": "locus_id",
        }
        for old, new in aliases.items():
            if old in frame and new not in frame:
                frame[new] = frame[old]
        identity = variant_column or (
            "variant_id.ss" if format == "susie" and "variant_id.ss" in frame else "variant_id"
        )
        if identity not in frame:
            raise ValueError(f"Missing variant column {identity!r}; specify variant_column")
        if effect_allele_column is not None:
            if effect_allele_column not in frame:
                raise ValueError(f"Missing effect_allele_column {effect_allele_column!r}")
            frame["effect_allele"] = frame[effect_allele_column]
        frame["input_variant_id"] = frame[identity]
        variants = []
        for value in frame[identity]:
            match = re.fullmatch(
                r"([^:]+):(\d+)\[b(37|38)\]([ACGTacgt]+),([ACGTacgt]+)", str(value)
            )
            if match is None:
                if not re.fullmatch(r"[^:]+:\d+:[ACGTacgt]+:[ACGTacgt]+", str(value)):
                    raise ValueError(
                        f"Unsupported Matos variant ID {value!r}; supply canonical coordinates"
                    )
                variants.append(str(value))
                continue
            chrom, pos, build, first, second = match.groups()
            if genome_build.lower() not in {f"grch{build}", "hg19" if build == "37" else "hg38"}:
                raise ValueError(
                    "Embedded variant build conflicts with genome_build; use an externally harmonized canonical file"
                )
            if allele_order not in {"ref-alt", "alt-ref"}:
                raise ValueError(
                    "Custom Matos IDs require explicit allele_order verified against the release VCF"
                )
            ref, alt = (first, second) if allele_order == "ref-alt" else (second, first)
            variants.append(f"{chrom}:{pos}:{ref}:{alt}")
        frame["variant_id"] = variants
    if context is not None:
        if "context" not in frame:
            frame["context"] = context
        else:
            frame["context"] = frame["context"].fillna(context)
    return normalize_regulatory_variants(
        frame, layer=layer, source=source, release=release, genome_build=genome_build
    )


def run_manifest(manifest_path: Path, outdir: Path) -> None:
    """Read local files, record byte checksums, and export deterministic result tables."""
    manifest = json.loads(manifest_path.read_text())
    if "rfu_qtl" not in manifest or set(manifest) - set(_TARGETS):
        raise ValueError("Manifest requires rfu_qtl and may contain eqtl, caqtl, gwas")
    tables, metadata = {}, {}
    for layer, spec in manifest.items():
        for field in ("source", "release", "genome_build", "files"):
            if not spec.get(field):
                raise ValueError(f"{layer} requires explicit {field}")
        frames, resources = [], []
        for entry in spec["files"]:
            path = (manifest_path.parent / entry["path"]).resolve()
            if not path.is_file():
                raise FileNotFoundError(path)
            digest = hashlib.sha256()
            with path.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
            sha256 = digest.hexdigest()
            if entry.get("sha256") and entry["sha256"] != sha256:
                raise ValueError(f"SHA256 mismatch for {path}")
            frame = pd.read_csv(path, sep=entry.get("sep", "\t"), dtype=object)
            adapted = adapt_matos_table(
                frame,
                layer=layer,
                **{k: spec[k] for k in ("source", "release", "genome_build")},
                context=entry.get("context", spec.get("context")),
                **{
                    k: entry[k]
                    for k in (
                        "format",
                        "allele_order",
                        "variant_column",
                        "effect_allele_column",
                        "column_map",
                    )
                    if k in entry
                },
            )
            adapted["sha256"] = sha256
            adapted["url"] = entry.get("url", spec.get("url"))
            frames.append(adapted)
            resources.append({**entry, "path": str(path), "sha256": sha256, "rows": len(frame)})
        tables[layer] = pd.concat(frames, ignore_index=True)
        metadata[layer] = {k: v for k, v in spec.items() if k != "files"}
        metadata[layer]["files"] = resources
    result = regulatory_triangulation(**tables, input_metadata=metadata)
    outputs = {
        "harmonized_rfu_qtl.tsv": result.harmonized_rfu_qtl,
        "regulatory_hits.tsv": result.matched_evidence,
        "rfu_regulatory_summary.tsv": result.rfu_summary,
        "variant_regulatory_summary.tsv": result.variant_summary,
        "unmatched_variants.tsv": result.unmatched_variants,
        "credible_set_overlap.tsv": result.credible_set_overlaps,
    }
    for name in [*outputs, "provenance.json"]:
        if (outdir / name).exists():
            raise FileExistsError(f"Refusing to overwrite {outdir / name}")
    outdir.mkdir(parents=True, exist_ok=True)
    for name, frame in outputs.items():
        frame.to_csv(outdir / name, sep="\t", index=False)
    (outdir / "provenance.json").write_text(
        json.dumps(result.provenance, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    run_manifest(args.manifest.resolve(), args.outdir)


if __name__ == "__main__":
    main()
