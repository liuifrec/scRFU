"""Offline, descriptive RFU genetic evidence integration; no causal inference."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from typing import Any

import numpy as np
import pandas as pd

from ._version import __version__

_TARGETS = {"rfu_qtl": "rfu_label", "eqtl": "gene", "caqtl": "peak", "gwas": "trait"}
_COORDS = ("chromosome", "position", "ref", "alt")
_OPTIONAL = (
    "variant_id",
    "rsid",
    "beta",
    "se",
    "pvalue",
    "allele_frequency",
    "pip",
    "odds_ratio",
    "effect_allele",
    "strand",
    "context",
    "credible_set_id",
    "locus_id",
    "release",
    "url",
    "sha256",
)
_FLAGS = ("has_rfu_qtl", "has_eqtl", "has_caqtl", "has_gwas", "eqtl_finemapped", "caqtl_finemapped")


@dataclass(frozen=True)
class RegulatoryEvidenceSchema:
    """Canonical columns for one evidence layer; coordinates may come from variant_id."""

    layer: str
    required: tuple[str, ...]
    optional: tuple[str, ...]


def regulatory_evidence_schema(layer: str) -> RegulatoryEvidenceSchema:
    """Return a schema for rfu_qtl, eqtl, caqtl, or gwas (no external resources)."""
    if layer not in _TARGETS:
        raise ValueError(f"Unknown evidence layer {layer!r}; choose from {tuple(_TARGETS)}")
    return RegulatoryEvidenceSchema(
        layer, (*_COORDS, "genome_build", _TARGETS[layer], "source"), _OPTIONAL
    )


@dataclass(frozen=True)
class RegulatoryTriangulationResult:
    """Pairwise evidence, association summaries, QC, set overlaps, and provenance."""

    harmonized_rfu_qtl: pd.DataFrame
    matched_evidence: pd.DataFrame
    unmatched_variants: pd.DataFrame
    rfu_summary: pd.DataFrame
    variant_summary: pd.DataFrame
    credible_set_overlaps: pd.DataFrame
    provenance: dict[str, Any]


def _text(value: Any) -> str | None:
    if pd.isna(value):
        return None
    value = str(value).strip()
    return value or None


def _build(value: Any) -> str | None:
    value = _text(value)
    aliases = {"hg19": "GRCh37", "grch37": "GRCh37", "hg38": "GRCh38", "grch38": "GRCh38"}
    return aliases.get(value.lower(), value) if value else None


def _chromosome(value: Any) -> str | None:
    value = _text(value)
    if value is None:
        return None
    value = re.sub(r"^chr", "", value, flags=re.IGNORECASE)
    if re.fullmatch(r"\d+(?:\.0+)?", value):
        return str(int(float(value)))
    if value.upper() in {"X", "Y", "M", "MT"}:
        return "MT" if value.upper() in {"M", "MT"} else value.upper()
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", value):
        raise ValueError(f"Invalid chromosome/contig: {value!r}")
    return value


def _signature(row: pd.Series) -> str:
    values = {
        str(k): None if pd.isna(v) else v.item() if isinstance(v, np.generic) else v
        for k, v in row.items()
    }
    return json.dumps(values, sort_keys=True, ensure_ascii=True, default=str)


def _ordered(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame.reset_index(drop=True)
    order = frame.apply(_signature, axis=1).sort_values(kind="stable").index
    return frame.loc[order].reset_index(drop=True)


def normalize_regulatory_variants(
    table: pd.DataFrame,
    *,
    layer: str,
    genome_build: str | None = None,
    source: str | None = None,
    release: str | None = None,
) -> pd.DataFrame:
    """Validate and normalize a copy of a canonical evidence table.

    Positions are one-based integers; alleles are forward-genomic A/C/G/T sequences.
    Coordinate IDs must be CHROM:POS:REF:ALT. rsIDs alone remain unresolved, never
    matched. Explicit effect_allele is required for direction comparisons; ALT is
    not assumed. No liftover, strand complementation, reference lookup, or indel
    left alignment is performed. Empty tables are accepted without columns.
    """
    schema = regulatory_evidence_schema(layer)
    if not isinstance(table, pd.DataFrame):
        raise TypeError("Evidence must be a pandas DataFrame")
    if not table.columns.is_unique:
        raise ValueError("Duplicate column names are not supported")
    frame = table.copy().reset_index(drop=True)
    for column, default in (
        ("genome_build", genome_build),
        ("source", source),
        ("release", release),
    ):
        if column not in frame:
            frame[column] = default
        elif default is not None:
            frame[column] = frame[column].fillna(default)
        if column == "genome_build" and default is not None:
            existing = frame[column].map(_build).dropna()
            if not existing.eq(_build(default)).all():
                raise ValueError("genome_build argument conflicts with table genome builds")
    target = _TARGETS[layer]
    if not frame.empty and target not in frame:
        raise ValueError(f"{layer} requires column {target!r}")
    if not frame.empty and not (
        set(_COORDS).issubset(frame) or "variant_id" in frame or "rsid" in frame
    ):
        raise ValueError(
            "Variant identity requires chromosome, position, ref, alt or variant_id/rsid"
        )
    for column in (*schema.required, *schema.optional):
        if column not in frame:
            frame[column] = pd.NA
    for column in (
        target,
        "source",
        "release",
        "url",
        "sha256",
        "context",
        "credible_set_id",
        "locus_id",
        "variant_id",
        "rsid",
    ):
        frame[column] = frame[column].map(_text)
    frame["genome_build"] = frame["genome_build"].map(_build)
    for column in (target, "source", "genome_build"):
        if frame[column].isna().any():
            raise ValueError(f"{layer} requires nonempty {column} for every record")
    for column in ("ref", "alt", "effect_allele"):
        frame[column] = frame[column].map(lambda x: s.upper() if (s := _text(x)) else None)
    frame["strand"] = frame["strand"].map(_text)
    if not frame["strand"].dropna().isin(["+", "-"]).all():
        raise ValueError("strand must be '+' or '-' when supplied")
    frame["chromosome"] = frame["chromosome"].map(_chromosome)
    # ID-based filling must also work when an input position column is nullable numeric.
    frame["position"] = frame["position"].astype(object)
    keys, positions, statuses = [], [], []
    for index, row in frame.iterrows():
        vid = row["variant_id"]
        parts = vid.split(":") if vid else []
        if len(parts) == 4:
            parsed = [_chromosome(parts[0]), parts[1], parts[2].upper(), parts[3].upper()]
            for col, val in zip(_COORDS, parsed, strict=True):
                if pd.isna(row[col]):
                    row[col] = val
                    frame.at[index, col] = val
                else:
                    same = str(row[col]) == val
                    if col == "position":
                        try:
                            same = float(row[col]) == float(val)
                        except (ValueError, TypeError):
                            same = False
                    if not same:
                        raise ValueError(f"variant_id conflicts with {col} at row {index}")
        if vid and re.fullmatch(r"rs\d+", vid, re.IGNORECASE):
            rsid = vid.lower()
            if row["rsid"] is not None and row["rsid"].lower() != rsid:
                raise ValueError("Conflicting rsid and variant_id")
            frame.at[index, "rsid"] = rsid
        if row["rsid"] is not None:
            if not re.fullmatch(r"rs\d+", row["rsid"], re.IGNORECASE):
                raise ValueError("rsid must have the form rs123")
            frame.at[index, "rsid"] = row["rsid"].lower()
        pos = row["position"]
        if not pd.isna(pos):
            try:
                number = Decimal(str(pos))
            except (InvalidOperation, ValueError, TypeError) as exc:
                raise ValueError("position must be a positive integer") from exc
            if (
                isinstance(pos, (bool, np.bool_))
                or not number.is_finite()
                or number < 1
                or number % 1
                or number > 2**53
            ):
                raise ValueError("position must be a positive integer <= 2**53")
            pos = int(number)
        positions.append(pos)
        for col in ("ref", "alt", "effect_allele"):
            if row[col] is not None and not re.fullmatch(r"[ACGT]+", row[col]):
                raise ValueError(
                    f"{col} must be an A/C/G/T sequence (split multiallelic records first)"
                )
        complete = not any(pd.isna(row[col]) for col in _COORDS)
        if complete and row["ref"] == row["alt"]:
            raise ValueError("ref and alt must differ")
        if (
            complete
            and row["effect_allele"] is not None
            and row["effect_allele"] not in (row["ref"], row["alt"])
        ):
            raise ValueError("effect_allele must equal ref or alt in the supplied orientation")
        keys.append(f"{row['chromosome']}:{pos}:{row['ref']}:{row['alt']}" if complete else None)
        statuses.append("resolved" if complete else "unresolved_variant_identity")
    frame["position"] = pd.array(positions, dtype="Int64")
    frame["variant_key"] = pd.array(keys, dtype="string")
    frame["identity_status"] = pd.array(statuses, dtype="string")
    for column in ("beta", "se", "pvalue", "allele_frequency", "pip", "odds_ratio"):
        values = pd.to_numeric(frame[column], errors="coerce").astype(float)
        if (frame[column].notna() & values.isna()).any() or np.isinf(values).any():
            raise ValueError(f"{column} must contain finite numbers or missing values")
        if (
            column in {"pvalue", "allele_frequency", "pip"}
            and not values.dropna().between(0, 1).all()
        ):
            raise ValueError(f"{column} must lie in [0, 1]")
        if column in {"se", "odds_ratio"} and (values.dropna() <= 0).any():
            raise ValueError(f"{column} must be positive")
        frame[column] = values
    if (frame["beta"].notna() & frame["odds_ratio"].notna()).any():
        raise ValueError("Supply beta or odds_ratio per row, not both; beta scales may differ")
    if layer != "gwas" and frame["odds_ratio"].notna().any():
        raise ValueError("odds_ratio is only supported for gwas")
    return _ordered(frame)


def _compatible(tables: Mapping[str, pd.DataFrame]) -> None:
    builds = {b for table in tables.values() for b in table["genome_build"].dropna()}
    if len(builds) > 1:
        raise ValueError(
            f"Incompatible genome builds: {sorted(builds)}; supply externally harmonized coordinates"
        )


def _threshold(value: float) -> None:
    if not np.isfinite(value) or not 0 <= value <= 1:
        raise ValueError("pip_threshold must lie in [0, 1]")


def _deduplicate(frame: pd.DataFrame, layer: str) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    frame = frame.copy()
    frame["record_id"] = [
        hashlib.sha256(_signature(row).encode()).hexdigest() for _, row in frame.iterrows()
    ]
    counts = frame["record_id"].value_counts()
    frame = frame.drop_duplicates("record_id").reset_index(drop=True)
    qc = [
        _qc(row, layer, "duplicate_record", int(counts[row["record_id"]] - 1))
        for _, row in frame.iterrows()
        if counts[row["record_id"]] > 1
    ]
    return frame, qc


def _qc(row: pd.Series, layer: str, reason: str, count: int = 1) -> dict[str, Any]:
    return {
        "layer": layer,
        "record_id": row["record_id"],
        "variant_id": row["variant_id"],
        "variant_key": row["variant_key"],
        "genome_build": row["genome_build"],
        "target": row[_TARGETS[layer]],
        "source": row["source"],
        "release": row["release"],
        "reason": reason,
        "record_count": count,
    }


def _direction(left: pd.Series, right: pd.Series) -> tuple[bool, Any, float, str]:
    if left["effect_allele"] is None or right["effect_allele"] is None:
        return False, pd.NA, np.nan, "missing_effect_allele"
    if left["strand"] == "-" or right["strand"] == "-":
        return False, pd.NA, np.nan, "unsupported_negative_strand"
    if {left["ref"], left["alt"]} in ({"A", "T"}, {"C", "G"}) and not (
        left["strand"] == right["strand"] == "+"
    ):
        return False, pd.NA, np.nan, "palindromic_strand_unknown"
    beta = right["beta"]
    if pd.isna(beta) and pd.notna(right["odds_ratio"]):
        beta = float(np.log(right["odds_ratio"]))
    aligned = beta if left["effect_allele"] == right["effect_allele"] else -beta
    if pd.isna(left["beta"]) or pd.isna(aligned):
        return True, pd.NA, aligned, "missing_effect_size"
    if left["beta"] == 0 or aligned == 0:
        return True, pd.NA, aligned, "zero_effect"
    return True, bool(np.sign(left["beta"]) == np.sign(aligned)), aligned, "aligned"


def _tier(has_eqtl: bool, has_caqtl: bool, has_gwas: bool) -> int:
    return (
        4
        if has_eqtl and has_caqtl and has_gwas
        else 3
        if has_eqtl and has_caqtl
        else 2
        if has_eqtl or has_caqtl
        else 1
    )


def _labels(series: pd.Series) -> str:
    return json.dumps(sorted(set(series.dropna().astype(str))), ensure_ascii=True)


def _summarize(associations: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    columns = [
        *by,
        "n_rfu_qtl_records",
        "n_variants",
        *_FLAGS,
        "evidence_count",
        "max_evidence_tier",
        "genes",
        "peaks",
        "traits",
    ]
    rows = []
    for key, group in associations.groupby(by, sort=True, dropna=False, observed=True):
        key = key if isinstance(key, tuple) else (key,)
        row = dict(zip(by, key, strict=True))
        row.update(n_rfu_qtl_records=len(group), n_variants=group["variant_key"].nunique())
        row.update({flag: bool(group[flag].any()) for flag in _FLAGS})
        row["evidence_count"] = sum(row[f"has_{layer}"] for layer in _TARGETS)
        row["max_evidence_tier"] = int(group["evidence_tier"].max())
        for label in ("genes", "peaks", "traits"):
            row[label] = json.dumps(sorted({s for cell in group[label] for s in json.loads(cell)}))
        rows.append(row)
    return pd.DataFrame(rows, columns=columns)


def credible_set_overlap(
    left: pd.DataFrame,
    right: pd.DataFrame,
    *,
    left_layer: str = "rfu_qtl",
    right_layer: str = "eqtl",
    pip_threshold: float = 0.95,
) -> pd.DataFrame:
    """Describe exact-key intersections of supplied credible-set memberships.

    Set IDs are scoped by target, context, source, release, locus_id and build.
    Only overlapping set pairs are returned. Sizes refer to supplied members;
    missing PIPs are not high-PIP support. This is not statistical colocalization.
    """
    _threshold(pip_threshold)
    a = normalize_regulatory_variants(left, layer=left_layer)
    b = normalize_regulatory_variants(right, layer=right_layer)
    _compatible({"left": a, "right": b})
    return _set_overlap(a, b, left_layer, right_layer, pip_threshold)


def _set_overlap(
    a: pd.DataFrame, b: pd.DataFrame, left_layer: str, right_layer: str, threshold: float
) -> pd.DataFrame:
    scope = ["genome_build", "source", "release", "context", "locus_id", "credible_set_id"]
    columns = [
        "left_layer",
        "right_layer",
        "left_set",
        "right_set",
        "left_size",
        "right_size",
        "n_shared",
        "shared_variants",
        "n_shared_high_pip",
        "shared_high_pip_variants",
        "jaccard",
    ]

    def groups(table: pd.DataFrame, layer: str) -> dict[str, dict[str, float]]:
        out = {}
        valid = table.loc[table["credible_set_id"].notna() & table["variant_key"].notna()]
        for key, group in valid.groupby(
            [*scope, _TARGETS[layer]], dropna=False, sort=True, observed=True
        ):
            identifier = json.dumps(
                {
                    c: None if pd.isna(v) else v
                    for c, v in zip([*scope, _TARGETS[layer]], key, strict=True)
                },
                sort_keys=True,
            )
            out[identifier] = group.groupby("variant_key")["pip"].max().to_dict()
        return out

    ag, bg = groups(a, left_layer), groups(b, right_layer)
    inverted: dict[str, set[str]] = {}
    for name, members in bg.items():
        for variant in members:
            inverted.setdefault(variant, set()).add(name)
    rows = []
    for aname, amembers in ag.items():
        candidates = {name for v in amembers for name in inverted.get(v, ())}
        for bname in sorted(candidates):
            bmembers = bg[bname]
            shared = sorted(amembers.keys() & bmembers.keys())
            high = [v for v in shared if amembers[v] >= threshold and bmembers[v] >= threshold]
            rows.append(
                dict(
                    left_layer=left_layer,
                    right_layer=right_layer,
                    left_set=aname,
                    right_set=bname,
                    left_size=len(amembers),
                    right_size=len(bmembers),
                    n_shared=len(shared),
                    shared_variants=json.dumps(shared),
                    n_shared_high_pip=len(high),
                    shared_high_pip_variants=json.dumps(high),
                    jaccard=len(shared) / len(amembers.keys() | bmembers.keys()),
                )
            )
    return pd.DataFrame(rows, columns=columns)


def regulatory_triangulation(
    rfu_qtl: pd.DataFrame,
    *,
    eqtl: pd.DataFrame | None = None,
    caqtl: pd.DataFrame | None = None,
    gwas: pd.DataFrame | None = None,
    pip_threshold: float = 0.95,
    allow_allele_reversal: bool = True,
    input_metadata: Mapping[str, Mapping[str, Any]] | None = None,
) -> RegulatoryTriangulationResult:
    """Integrate user-supplied RFU-QTL, regulatory QTL and GWAS evidence offline.

    No significance filtering is implicit. Each matched row is an RFU association
    paired with one external evidence record (no cross-layer Cartesian products).
    Reversed allele pairs are flagged separately from exact keys. Summary tiers
    are descriptive; GWAS alone leaves an association at tier 1. All inputs must
    declare the same genome build. Metadata can supply source/release/build defaults.
    """
    _threshold(pip_threshold)
    inputs = {"rfu_qtl": rfu_qtl, "eqtl": eqtl, "caqtl": caqtl, "gwas": gwas}
    metadata = dict(input_metadata or {})
    if set(metadata) - set(inputs):
        raise ValueError("input_metadata keys must be evidence layer names")
    tables, qc, provenance_inputs = {}, [], {}
    for layer, table in inputs.items():
        meta = dict(metadata.get(layer, {}))
        frame = normalize_regulatory_variants(
            table if table is not None else pd.DataFrame(),
            layer=layer,
            **{k: meta[k] for k in ("genome_build", "source", "release") if k in meta},
        )
        frame, duplicates = _deduplicate(frame, layer)
        tables[layer] = frame
        qc.extend(duplicates)
        provenance_inputs[layer] = {
            "supplied": table is not None,
            "input_rows": 0 if table is None else len(table),
            "unique_rows": len(frame),
            "metadata": meta,
            "resources": frame[["source", "release", "url", "sha256", "genome_build"]]
            .drop_duplicates()
            .to_dict("records"),
        }
    _compatible(tables)
    base = tables["rfu_qtl"].copy()
    for flag in _FLAGS:
        base[flag] = flag == "has_rfu_qtl"
    for label in ("genes", "peaks", "traits"):
        base[label] = "[]"
    base["shared_variant_exact"] = False
    matched = []
    used: dict[str, set[str]] = {layer: set() for layer in inputs}
    indexes = {}
    for layer, table in tables.items():
        if layer == "rfu_qtl":
            continue
        indexes[layer] = {
            key: group
            for key, group in table.loc[table["variant_key"].notna()].groupby(
                "variant_key", sort=False
            )
        }
    for index, rfu in base.iterrows():
        if pd.isna(rfu["variant_key"]):
            continue
        reverse = f"{rfu['chromosome']}:{rfu['position']}:{rfu['alt']}:{rfu['ref']}"
        for layer in ("eqtl", "caqtl", "gwas"):
            hits = []
            for key in [rfu["variant_key"], *([reverse] if allow_allele_reversal else [])]:
                if key in indexes[layer]:
                    hits.extend(row for _, row in indexes[layer][key].iterrows())
            if not hits:
                continue
            base.at[index, f"has_{layer}"] = True
            base.at[index, {"eqtl": "genes", "caqtl": "peaks", "gwas": "traits"}[layer]] = _labels(
                pd.Series([h[_TARGETS[layer]] for h in hits])
            )
            for hit in hits:
                exact = rfu["variant_key"] == hit["variant_key"]
                harmonized, concordant, aligned, reason = _direction(rfu, hit)
                finemapped = bool(hit["pip"] >= pip_threshold)
                if layer != "gwas" and finemapped:
                    base.at[index, f"{layer}_finemapped"] = True
                if exact:
                    base.at[index, "shared_variant_exact"] = True
                used["rfu_qtl"].add(rfu["record_id"])
                used[layer].add(hit["record_id"])
                record = {
                    "rfu_record_id": rfu["record_id"],
                    "rfu_label": rfu["rfu_label"],
                    "variant_key": rfu["variant_key"],
                    "genome_build": rfu["genome_build"],
                    "evidence_layer": layer,
                    "shared_variant_exact": exact,
                    "alleles_reversed": not exact,
                    "allele_harmonized": harmonized,
                    "direction_concordant": concordant,
                    "direction_status": reason,
                    "aligned_evidence_beta": aligned,
                    "rfu_beta": rfu["beta"],
                    "rfu_effect_allele": rfu["effect_allele"],
                    "rfu_pip": rfu["pip"],
                    "evidence_finemapped": finemapped,
                    "shared_high_pip": bool(rfu["pip"] >= pip_threshold and finemapped),
                }
                record.update({f"evidence_{col}": value for col, value in hit.items()})
                matched.append(record)
    base["evidence_count"] = base[[f"has_{layer}" for layer in inputs]].sum(axis=1).astype(int)
    base["evidence_tier"] = [_tier(r.has_eqtl, r.has_caqtl, r.has_gwas) for r in base.itertuples()]
    for layer, table in tables.items():
        for _, row in table.iterrows():
            if row["record_id"] not in used[layer]:
                reason = (
                    "unresolved_variant_identity"
                    if pd.isna(row["variant_key"])
                    else "no_matching_evidence"
                    if layer == "rfu_qtl"
                    else "no_matching_rfu_qtl"
                )
                qc.append(_qc(row, layer, reason))
    match_columns = [
        "rfu_record_id",
        "rfu_label",
        "variant_key",
        "genome_build",
        "evidence_layer",
        "shared_variant_exact",
        "alleles_reversed",
        "allele_harmonized",
        "direction_concordant",
        "direction_status",
        "aligned_evidence_beta",
        "rfu_beta",
        "rfu_effect_allele",
        "rfu_pip",
        "evidence_finemapped",
        "shared_high_pip",
    ]
    evidence_columns = sorted(
        {f"evidence_{c}" for layer in ("eqtl", "caqtl", "gwas") for c in tables[layer].columns}
    )
    matches = _ordered(pd.DataFrame(matched, columns=[*match_columns, *evidence_columns]))
    for col in (
        "shared_variant_exact",
        "alleles_reversed",
        "allele_harmonized",
        "direction_concordant",
        "evidence_finemapped",
        "shared_high_pip",
    ):
        matches[col] = matches[col].astype("boolean")
    qc_table = _ordered(
        pd.DataFrame(
            qc,
            columns=[
                "layer",
                "record_id",
                "variant_id",
                "variant_key",
                "genome_build",
                "target",
                "source",
                "release",
                "reason",
                "record_count",
            ],
        )
    )
    overlap_tables = [
        _set_overlap(tables["rfu_qtl"], tables[layer], "rfu_qtl", layer, pip_threshold)
        for layer in ("eqtl", "caqtl", "gwas")
    ]
    nonempty_overlaps = [table for table in overlap_tables if not table.empty]
    overlaps = (
        pd.concat(nonempty_overlaps, ignore_index=True) if nonempty_overlaps else overlap_tables[0]
    )
    provenance = {
        "schema_version": 1,
        "scrfu_version": __version__,
        "method": "descriptive_regulatory_triangulation",
        "parameters": {
            "pip_threshold": pip_threshold,
            "allow_allele_reversal": allow_allele_reversal,
        },
        "inputs": provenance_inputs,
        "significance_filter": None,
        "formal_colocalization": False,
        "causal_inference": False,
        "duplicate_policy": "collapse identical normalized records; retain distinct associations",
        "variant_policy": "build-aware exact keys or explicit reversed allele pairs; no liftover or strand flips",
    }
    return RegulatoryTriangulationResult(
        base,
        matches,
        qc_table,
        _summarize(base, ["rfu_label"]),
        _summarize(base.loc[base["variant_key"].notna()], ["genome_build", "variant_key"]),
        overlaps,
        provenance,
    )


def join_regulatory_summary(
    summary: pd.DataFrame, result: RegulatoryTriangulationResult
) -> pd.DataFrame:
    """Left join RFU-level evidence on rfu_label, preserving input order and index.

    Repeated RFU labels (e.g. one per donor) are allowed on the left. Unobserved
    RFUs retain missing evidence, not false evidence. Column collisions raise.
    """
    if "rfu_label" not in summary or summary["rfu_label"].isna().any():
        raise ValueError("summary requires nonmissing rfu_label")
    collision = (set(summary) & set(result.rfu_summary)) - {"rfu_label"}
    if collision:
        raise ValueError(f"Regulatory summary column collision: {sorted(collision)}")
    left = summary.copy()
    left["rfu_label"] = left["rfu_label"].map(_text)
    if left["rfu_label"].isna().any():
        raise ValueError("summary requires nonempty rfu_label")
    joined = left.merge(
        result.rfu_summary, on="rfu_label", how="left", sort=False, validate="many_to_one"
    )
    joined.index = summary.index.copy()
    joined["rfu_label"] = summary["rfu_label"].to_numpy()
    return joined
