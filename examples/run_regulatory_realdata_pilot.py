"""Analyze the prepared RfuWAS × Matos GRCh38 exact-variant lookup offline.

Uses a single Bonferroni family across all retrieved variant–target tests in
both molecular layers, plus separately identified released independent QTLs.
This analyst-defined filter does not reconstruct unreleased source FDR values.
Run matos_regulatory_lookup.py for each layer first. Outputs stay outside Git.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from scrfu.tl import regulatory_triangulation


def labels(values: pd.Series) -> str:
    return json.dumps(sorted(set(values.dropna().astype(str))))


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep="\t",
        dtype={
            "rfu_label": str,
            "chromosome": str,
            "phecode": str,
            "strongest_disease_phecode": str,
            "autoimmune_disease_phecode": str,
        },
    )


def save(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False)


def evidence_sets(raw: pd.DataFrame, independent: pd.DataFrame, threshold: float) -> dict:
    return {
        "tested": set(raw.variant_key),
        "nominal": set(raw.loc[raw.pvalue < 0.05, "variant_key"]),
        "lookup_bonferroni": set(raw.loc[raw.pvalue < threshold, "variant_key"]),
        "independent": set(independent.variant_key),
        "evidence": set(raw.loc[raw.pvalue < threshold, "variant_key"])
        | set(independent.variant_key),
    }


def run(root: Path, *, partial_caqtl: bool = False) -> None:
    prepared, baseline = root / "prepared", root / "results"
    out = baseline / "partial_caqtl_checkpoint" if partial_caqtl else baseline
    out.mkdir(parents=True, exist_ok=True)
    rfu = read(prepared / "rfuwas_data1_rfuqtl_grch38.tsv")
    raw = {
        layer: read(
            prepared
            / (
                "matos_caqtl_chr7_partial_at_rfuqtl_variants.tsv"
                if partial_caqtl and layer == "caqtl"
                else f"matos_{layer}_at_rfuqtl_variants.tsv"
            )
        )
        for layer in ("eqtl", "caqtl")
    }
    all_independent = {
        layer: read(
            prepared
            / (
                "matos_caqtl_chr6_partial_independent.tsv"
                if partial_caqtl and layer == "caqtl"
                else f"matos_{layer}_independent_chr6_chr7.tsv"
            )
        )
        for layer in raw
    }
    if partial_caqtl and "source_file" not in all_independent["caqtl"]:
        all_independent["caqtl"]["source_file"] = (
            "cd4_qsmooth_cpm_chromatin_narrowpeaks.independent_cis_QTL_pairs.chr6.csv"
        )
    independent = {
        layer: frame.loc[frame.variant_key.isin(rfu.variant_key)].copy()
        for layer, frame in all_independent.items()
    }
    family_size = sum(len(frame) for frame in raw.values())
    threshold = 0.05 / family_size if family_size else 0.0
    counts = {
        "rfu_qtl_associations": len(rfu),
        "rfu_qtl_variants": rfu.variant_key.nunique(),
        "rfus": rfu.rfu_label.nunique(),
        "lookup_test_family_size": family_size,
        "lookup_bonferroni_threshold": threshold,
        "source_marginal_significance": "unknown: permutation/FDR tables absent from archives",
        "partial_caqtl": partial_caqtl,
        "caqtl_scope": "chr7 nominal and chr6 independent only; archive checksum unverified"
        if partial_caqtl
        else "chr6 and chr7; verified archive",
    }
    tables, sets = {}, {}
    for layer, frame in raw.items():
        target = "gene" if layer == "eqtl" else "peak"
        if frame.duplicated(["variant_key", target]).any():
            raise ValueError(
                "Duplicate variant–target tests; establish test family before correction"
            )
        frame["lookup_bonferroni"] = frame.pvalue < threshold
        frame["lookup_bonferroni_threshold"] = threshold
        indep = independent[layer]
        conditional_cols = ["variant_key", target, "pvalue", "beta", "se", "rank", "source_file"]
        frame = frame.merge(
            indep[conditional_cols].rename(
                columns={c: f"conditional_{c}" for c in conditional_cols[2:]}
            ),
            on=["variant_key", target],
            how="left",
            validate="one_to_one",
        )
        selected = frame.loc[
            frame.lookup_bonferroni | frame.independent.astype("boolean").fillna(False)
        ].copy()
        selected["analysis_evidence_basis"] = "lookup_bonferroni_or_released_independent"
        covered = set(zip(frame.variant_key, frame[target], strict=True))
        additional = indep.loc[
            [(v, t) not in covered for v, t in zip(indep.variant_key, indep[target], strict=True)]
        ].copy()
        if not additional.empty:
            additional["analysis_evidence_basis"] = "released_independent_only"
            selected = pd.concat([selected, additional], ignore_index=True)
        tables[layer] = selected
        sets[layer] = evidence_sets(frame, indep, threshold)
        save(frame, out / f"matos_{layer}_lookup_classified.tsv")
        save(indep, out / f"matos_{layer}_independent_exact.tsv")
        counts[layer] = {
            "tested_variant_target_pairs": len(frame),
            "independent_variant_target_pairs": len(indep),
            "source_significant_total": None,
        }
        for category, variants in sets[layer].items():
            subset = rfu.loc[rfu.variant_key.isin(variants)]
            counts[layer][category] = {
                "variants": len(variants),
                "rfus": subset.rfu_label.nunique(),
                "rfu_qtl_associations": len(subset),
            }
    both = sets["eqtl"]["evidence"] & sets["caqtl"]["evidence"]
    counts["both_evidence"] = {
        "variants": len(both),
        "rfus": rfu.loc[rfu.variant_key.isin(both), "rfu_label"].nunique(),
    }
    both_independent = sets["eqtl"]["independent"] & sets["caqtl"]["independent"]
    counts["both_independent"] = {
        "variants": len(both_independent),
        "rfus": rfu.loc[rfu.variant_key.isin(both_independent), "rfu_label"].nunique(),
    }
    if partial_caqtl:
        counts["caqtl"]["independent_chr6_observed"] = counts["caqtl"]["independent"]
        counts["caqtl"]["independent_variant_target_pairs_chr6_observed"] = counts["caqtl"][
            "independent_variant_target_pairs"
        ]
        counts["caqtl"]["independent_variant_target_pairs"] = None
        counts["caqtl"]["independent"] = {
            "variants": None,
            "rfus": None,
            "rfu_qtl_associations": None,
        }
        counts["both_independent"] = {"variants": None, "rfus": None}
    metadata = {
        "analysis": "GRCh38 exact overlap; allele reversal disabled; effects unresolved",
        "evidence_filter": "lookup p < 0.05 / all retrieved eQTL+caQTL tests OR released independent QTL",
        "counts": counts,
        "archives": {
            layer: json.loads(
                (
                    prepared
                    / (
                        "matos_caqtl_partial_manifest.json"
                        if partial_caqtl and layer == "caqtl"
                        else f"matos_{layer}_archive_manifest.json"
                    )
                ).read_text()
            )
            for layer in raw
        },
        "rfuwas": json.loads((baseline / "rfuwas_context_provenance.json").read_text()),
    }
    print(json.dumps(counts, indent=2), flush=True)
    result = regulatory_triangulation(
        rfu,
        **tables,
        allow_allele_reversal=False,
        input_metadata={
            "rfu_qtl": metadata["rfuwas"],
            **{
                layer: {
                    "archive": metadata["archives"][layer],
                    "evidence_filter": metadata["evidence_filter"],
                }
                for layer in raw
            },
        },
    )
    assert result.matched_evidence.direction_concordant.isna().all()
    for name, frame in {
        "regulatory_hits.tsv": result.matched_evidence,
        "rfu_regulatory_summary.tsv": result.rfu_summary,
        "variant_regulatory_summary.tsv": result.variant_summary,
        "unmatched_variants.tsv": result.unmatched_variants,
        "harmonized_rfu_qtl.tsv": result.harmonized_rfu_qtl,
    }.items():
        save(frame, out / name)
    provenance = result.provenance
    provenance["pilot"] = metadata
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    candidates = result.harmonized_rfu_qtl.copy()
    for layer in raw:
        for category, variants in sets[layer].items():
            candidates[f"matos_{layer}_{category}"] = candidates.variant_key.isin(variants)
        target = "gene" if layer == "eqtl" else "peak"
        summary = (
            tables[layer]
            .groupby("variant_key", as_index=False)
            .agg(**{f"{layer}_min_pvalue": ("pvalue", "min"), f"{layer}_targets": (target, labels)})
        )
        candidates = candidates.merge(summary, on="variant_key", how="left", validate="many_to_one")
    context = read(baseline / "rfuqtl_cellstate_context.tsv")
    candidates = candidates.merge(context, on="rfu_label", how="left", validate="many_to_one")
    supported = candidates.loc[candidates.has_eqtl | candidates.has_caqtl].copy()
    save(
        result.matched_evidence.merge(context, on="rfu_label", how="left", validate="many_to_one"),
        out / "regulatory_hits_with_cellstate.tsv",
    )
    disease = read(prepared / "rfuwas_data4_annotations.tsv")
    disease = disease.loc[disease.rfu_label.isin(supported.rfu_label)].sort_values(
        ["pvalue", "rfu_label", "phecode"]
    )
    save(disease, out / "regulatory_rfu_disease_links.tsv")
    disease_summary = disease.groupby("rfu_label", as_index=False).agg(
        n_phecode_associations=("phecode", "nunique")
    )
    for prefix, subset in (
        ("strongest", disease),
        ("autoimmune", disease.loc[disease.group == "autoimmune"]),
    ):
        best = subset.drop_duplicates("rfu_label")[
            ["rfu_label", "description", "phecode", "pvalue", "effect"]
        ]
        best = best.rename(columns={c: f"{prefix}_disease_{c}" for c in best if c != "rfu_label"})
        disease_summary = disease_summary.merge(
            best, on="rfu_label", how="left", validate="one_to_one"
        )
    supported = supported.merge(disease_summary, on="rfu_label", how="left", validate="many_to_one")
    supported["has_cd4_annotation"] = supported.enrichment == "CD4-enriched"
    supported["has_cellstate_annotation"] = supported.Enriched_Group.notna()
    supported["has_disease_association"] = supported.n_phecode_associations.fillna(0) > 0
    support_flags = [
        "has_rfu_qtl",
        "has_eqtl",
        "has_caqtl",
        "matos_eqtl_independent",
        "matos_caqtl_independent",
        "has_cd4_annotation",
        "has_cellstate_annotation",
        "has_disease_association",
    ]
    supported["candidate_evidence_count"] = supported[support_flags].sum(axis=1)
    if partial_caqtl:
        supported["matos_caqtl_independent"] = pd.NA
        supported["candidate_evidence_count_is_lower_bound"] = True
    supported = supported.sort_values(
        ["candidate_evidence_count", "pvalue", "variant_key", "rfu_label"],
        ascending=[False, True, True, True],
    )
    save(supported, out / "ranked_regulatory_candidates.tsv")
    unique_rfus = supported.drop_duplicates("rfu_label")
    counts["supported_rfu_context"] = {
        "cd4_cd8": unique_rfus.enrichment.fillna("unannotated").value_counts().to_dict(),
        "cellstate": unique_rfus.Enriched_Group.fillna("unannotated").value_counts().to_dict(),
        "disease_rfus": disease.rfu_label.nunique(),
        "disease_links": len(disease),
    }
    # Proximity is only explored if fewer than ten exact supported variants exist.
    if len(sets["eqtl"]["evidence"] | sets["caqtl"]["evidence"]) < 10:
        proximity = []
        for layer, frame in all_independent.items():
            target = "gene" if layer == "eqtl" else "peak"
            for variant in rfu.drop_duplicates("variant_key").itertuples():
                near = frame.loc[
                    (frame.chromosome == variant.chromosome)
                    & (frame.position - variant.position).abs().le(250_000)
                ].copy()
                near["distance_bp"] = (near.position - variant.position).abs()
                near["within_100kb"] = near.distance_bp <= 100_000
                near["rfu_variant_key"] = variant.variant_key
                near["evidence_layer"] = layer
                near["target"] = near[target]
                near["match_type"] = "proximity_only"
                near = near.loc[near.variant_key != variant.variant_key]
                proximity.append(near)
        save(pd.concat(proximity, ignore_index=True), out / "proximity_candidates.tsv")
    else:
        counts["proximity_analysis"] = (
            "not triggered: at least ten variants have exact regulatory support"
        )
    (out / "overlap_counts.json").write_text(json.dumps(counts, indent=2) + "\n")
    plot_coverage(result.harmonized_rfu_qtl, out, provisional=partial_caqtl)
    print("Saved pilot outputs", flush=True)


def plot_coverage(associations: pd.DataFrame, out: Path, *, provisional: bool = False) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    variants = associations.drop_duplicates("variant_key")
    categories = {
        "None retrieved": (~variants.has_eqtl & ~variants.has_caqtl).sum(),
        "eQTL only": (variants.has_eqtl & ~variants.has_caqtl).sum(),
        "caQTL only": (~variants.has_eqtl & variants.has_caqtl).sum(),
        "Both": (variants.has_eqtl & variants.has_caqtl).sum(),
    }
    fig, ax = plt.subplots(figsize=(6.2, 3.8), layout="constrained")
    bars = ax.bar(
        categories.keys(), categories.values(), color=["#999999", "#4C78A8", "#F58518", "#54A24B"]
    )
    ax.bar_label(bars)
    ax.set(
        ylabel="Unique RFU-QTL variants",
        title=(
            "Provisional: caQTL chr7 only, checksum pending"
            if provisional
            else "Exact CD4 regulatory evidence at RFU-QTL variants"
        ),
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(out / "regulatory_evidence_coverage.png", dpi=200)
    fig.savefig(out / "regulatory_evidence_coverage.pdf")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument(
        "--partial-caqtl",
        action="store_true",
        help="Use recovered chr7 nominal / chr6 independent caQTL files; provisional outputs only",
    )
    args = parser.parse_args()
    run(args.root, partial_caqtl=args.partial_caqtl)
