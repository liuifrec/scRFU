"""Reproduce RFU-QTL input QC and RFU-level context before Matos intersection.

Requires local RfuWAS XLSX (openpyxl) and UCSC GRCh38 sequence JSON intervals.
These outputs describe all Data 1 RFUs, never a regulatory-supported subset.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

try:
    from .rfuwas_regulatory_prepare import prepare_rfuwas_data1
except ImportError:
    from rfuwas_regulatory_prepare import prepare_rfuwas_data1

WORKBOOK_SHA256 = "42742f4a30548c1184f8de022a2adc7b6b9abd321e81a01ea740fad5527891d9"


def reference_qc(variants: pd.DataFrame, references: list[dict]) -> pd.DataFrame:
    """Compare every distinct variant REF with covering GRCh38 sequence; never fix it."""
    for reference in references:
        if reference["genome"] != "hg38":
            raise ValueError("Expected hg38 reference sequences")
        if len(reference["dna"]) != reference["end"] - reference["start"]:
            raise ValueError("Reference sequence length disagrees with coordinates")
    rows = []
    for row in variants.drop_duplicates("variant_key").itertuples():
        if row.genome_build != "GRCh38":
            raise ValueError("Expected GRCh38 variants")
        start = int(row.position) - 1
        end = start + len(row.ref)
        bases = {
            r["dna"][start - r["start"] : end - r["start"]].upper()
            for r in references
            if r["chrom"].removeprefix("chr") == str(row.chromosome)
            and r["start"] <= start
            and r["end"] >= end
        }
        if len(bases) > 1:
            raise ValueError("Conflicting reference sequences")
        sequence = next(iter(bases), None)
        rows.append(
            {
                "variant_key": row.variant_key,
                "genome_build": row.genome_build,
                "ref": row.ref,
                "reference_sequence": sequence,
                "reference_checked": sequence is not None,
                "reference_match": sequence == row.ref if sequence is not None else pd.NA,
            }
        )
    return pd.DataFrame(rows)


def run(workbook: Path, root: Path, references: list[Path]) -> None:
    """Prepare Data 1–4, allele QC, baseline annotation tables, and provenance."""
    digest = hashlib.sha256(workbook.read_bytes()).hexdigest()
    if digest != WORKBOOK_SHA256:
        raise ValueError("Workbook SHA256 differs from the audited publication release")
    prepared, results = root / "prepared", root / "results"
    prepared.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)
    rfu = prepare_rfuwas_data1(
        pd.read_excel(workbook, sheet_name="Data 1", header=1), release="s42003-024-07010-x"
    )
    rfu["sha256"] = digest
    rfu.to_csv(prepared / "rfuwas_data1_rfuqtl_grch38.tsv", sep="\t", index=False)
    annotations = {}
    for number in (2, 3, 4):
        table = pd.read_excel(workbook, sheet_name=f"Data {number}", header=1)
        table["rfu_label"] = table["RFU"].map(lambda value: str(int(value)))
        table.to_csv(prepared / f"rfuwas_data{number}_annotations.tsv", sep="\t", index=False)
        annotations[number] = table
    ref = reference_qc(rfu, [json.loads(path.read_text()) for path in references])
    ref.to_csv(results / "rfuwas_reference_qc.tsv", sep="\t", index=False)
    context = rfu.groupby("rfu_label", as_index=False).agg(
        n_rfuqtl_variants=("variant_key", "nunique"), min_rfuqtl_pvalue=("pvalue", "min")
    )
    for number in (2, 3):
        table = annotations[number].drop(columns="RFU")
        context = context.merge(table, on="rfu_label", how="left", validate="one_to_one")
    context.to_csv(results / "rfuqtl_cellstate_context.tsv", sep="\t", index=False)
    disease = annotations[4].loc[annotations[4].rfu_label.isin(rfu.rfu_label)]
    disease.to_csv(results / "rfuqtl_disease_context.tsv", sep="\t", index=False)
    qc = {
        "association_count": len(rfu),
        "unique_variants": rfu.variant_key.nunique(),
        "unique_rfus": rfu.rfu_label.nunique(),
        "duplicate_variant_rfu_pairs": int(rfu.duplicated(["variant_key", "rfu_label"]).sum()),
        "unresolved_variants": int(rfu.variant_key.isna().sum()),
        "chromosome_associations": rfu.chromosome.value_counts().to_dict(),
        "trb_locus_interval": "GRCh38 chr7:142299011-142813287 (1-based inclusive)",
        "trb_locus_associations": int(
            ((rfu.chromosome == "7") & rfu.position.between(142299011, 142813287)).sum()
        ),
        "trb_flanking_interval": "GRCh38 chr7:141299011-143813287 (1-based inclusive)",
        "trb_flanking_associations": int(
            ((rfu.chromosome == "7") & rfu.position.between(141299011, 143813287)).sum()
        ),
        "reference_n_checked": int(ref.reference_checked.sum()),
        "reference_n_match": int(ref.reference_match.fillna(False).sum()),
        "reference_n_mismatch": int((ref.reference_match == False).sum()),  # noqa: E712
        "baseline_cd4_cd8": context.enrichment.fillna("unannotated").value_counts().to_dict(),
        "baseline_cellstate": context.Enriched_Group.fillna("unannotated").value_counts().to_dict(),
        "baseline_disease_links": len(disease),
        "baseline_disease_rfus": disease.rfu_label.nunique(),
    }
    qc["reference_match_fraction"] = (
        qc["reference_n_match"] / qc["reference_n_checked"] if qc["reference_n_checked"] else None
    )
    (results / "rfuwas_qc.json").write_text(json.dumps(qc, indent=2) + "\n")
    provenance = {
        "workbook": str(workbook.resolve()),
        "sha256": digest,
        "references": [
            {"path": str(p.resolve()), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
            for p in references
        ],
        "genome_build": "GRCh38",
        "effect_orientation": "unresolved",
        "scope": "all RFU-QTL RFUs before Matos filtering; not regulatory-supported RFUs",
    }
    (results / "rfuwas_context_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", required=True, type=Path)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--reference-json", required=True, nargs="+", type=Path)
    args = parser.parse_args()
    run(args.workbook, args.root, args.reference_json)
