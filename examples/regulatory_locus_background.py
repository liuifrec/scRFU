"""Exploratory Matos-only TRB background, without enrichment p-values.

Stream verified chr7 nominal members, summarize per variant, and standardize
control summaries to query strata. This is conditional on Matos coverage and
does not supply the unobserved RfuWAS tested universe or account for LD.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

try:
    from .matos_regulatory_lookup import digest
    from .regulatory_manuscript_audit import TRB
    from .run_regulatory_realdata_pilot import read, save
except ImportError:
    from matos_regulatory_lookup import digest
    from regulatory_manuscript_audit import TRB
    from run_regulatory_realdata_pilot import read, save


def summarize(path: Path, independent: pd.DataFrame) -> pd.DataFrame:
    retained = []
    for batch in pq.ParquetFile(path).iter_batches(
        batch_size=100_000,
        columns=["variant_id", "phenotype_id", "af", "pval_nominal", "start_distance"],
    ):
        frame = batch.to_pandas()
        frame["position"] = frame.variant_id.str.extract(r"^(?:chr)?7:(\d+)", expand=False).astype(
            int
        )
        frame = frame.loc[frame.position.between(TRB[0] - 250_000, TRB[1] + 250_000)]
        if not frame.empty:
            retained.append(frame)
    frame = pd.concat(retained, ignore_index=True)
    if frame.duplicated(["variant_id", "phenotype_id"]).any():
        raise ValueError("Repeated variant-target rows would inflate target opportunity")
    frame["abs_distance"] = frame.start_distance.abs()
    summary = frame.groupby("variant_id", as_index=False).agg(
        position=("position", "first"),
        af=("af", "first"),
        min_p=("pval_nominal", "min"),
        n_targets=("phenotype_id", "nunique"),
        min_target_distance=("abs_distance", "min"),
    )
    alleles = summary.variant_id.str.extract(r"^(?:chr)?7:\d+\[b38\]([ACGT]+),([ACGT]+)$")
    if alleles.isna().any().any():
        raise ValueError("Unexpected source identity")
    summary["variant_key"] = (
        "7:" + summary.position.astype(str) + ":" + alleles[0] + ":" + alleles[1]
    )
    summary["variant_class"] = np.where(
        (alleles[0].str.len() == 1) & (alleles[1].str.len() == 1), "SNV", "indel"
    )
    summary["maf"] = np.minimum(summary.af, 1 - summary.af)
    summary["independent"] = summary.variant_key.isin(independent.variant_key)
    summary["min_logp"] = -np.log10(summary.min_p.clip(lower=np.finfo(float).tiny))
    return summary


def run(root: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    rfu = read(root / "prepared/rfuwas_data1_rfuqtl_grch38.tsv")
    keys = set(rfu.variant_key)
    threshold = json.loads((out / "overlap_counts.json").read_text())["lookup_bonferroni_threshold"]
    provenance = {
        "scope": "exploratory_matos_locus_background",
        "TRB": TRB,
        "flank_bp": 250_000,
        "position_bin_bp": 100_000,
        "maf_bin_width": 0.1,
        "minimum_controls_per_query": 5,
        "pvalues": "none: LD/exchangeability/RfuWAS universe unresolved",
    }
    layers, summaries = {}, []
    for layer in ("eqtl", "caqtl"):
        manifest = json.loads((root / f"prepared/matos_{layer}_archive_manifest.json").read_text())
        entry = next(x for x in manifest["selected"] if x["filename"].endswith("chr7.parquet"))
        path = root / f"prepared/matos_{layer}_selected" / entry["filename"]
        if digest(path) != entry["sha256"]:
            raise ValueError("Verified source member changed")
        independent = read(root / f"prepared/matos_{layer}_independent_chr6_chr7.tsv")
        cache = out / f"{layer}_trb_variant_background.tsv"
        fingerprint = {
            "nominal_sha256": entry["sha256"],
            "independent_sha256": digest(
                root / f"prepared/matos_{layer}_independent_chr6_chr7.tsv"
            ),
            "settings": provenance,
            "script_sha256": digest(Path(__file__)),
        }
        cache_meta = cache.with_suffix(".json")
        if cache_meta.exists() and json.loads(cache_meta.read_text()) == fingerprint:
            frame = read(cache)
        else:
            frame = summarize(path, independent)
            save(frame, cache)
            cache_meta.write_text(json.dumps(fingerprint, indent=2) + "\n")
        frame["query"] = frame.variant_key.isin(keys)
        frame["query_threshold_pass"] = frame.min_p < threshold
        frame["position_bin"] = frame.position // 100_000
        frame["maf_bin"] = (frame.maf / 0.1).astype(int)
        frame["targets_bin"] = pd.cut(
            frame.n_targets, [0, 10, 25, 50, 100, 200, 500, 10000]
        ).astype(str)
        frame["distance_bin"] = pd.cut(
            frame.min_target_distance, [-1, 1000, 10000, 100000, 250000, 1000000]
        ).astype(str)
        strata = ["variant_class", "position_bin", "maf_bin", "targets_bin", "distance_bin"]
        controls = (
            frame.loc[~frame["query"]]
            .groupby(strata)
            .agg(
                controls=("variant_key", "size"),
                control_min_logp=("min_logp", "mean"),
                control_independent_rate=("independent", "mean"),
                control_query_threshold_rate=("query_threshold_pass", "mean"),
            )
        )
        matched = frame.loc[frame["query"]].merge(
            controls, on=strata, how="left", validate="many_to_one"
        )
        matched["usable_background"] = matched.controls.ge(5)
        matched["evidence_scope"] = "exploratory_matos_locus_background"
        save(matched, out / f"{layer}_matched_locus_background.tsv")
        usable = matched.loc[matched.usable_background]
        summaries.append(
            {
                "layer": layer,
                "query_variants_in_window": int(frame["query"].sum()),
                "control_variants_in_window": int((~frame["query"]).sum()),
                "matched_query_variants": len(usable),
                "observed_mean_min_logp": usable.min_logp.mean(),
                "control_standardized_mean_min_logp": usable.control_min_logp.mean(),
                "observed_independent_rate": usable.independent.mean(),
                "control_standardized_independent_rate": usable.control_independent_rate.mean(),
                "observed_query_threshold_rate": usable.query_threshold_pass.mean(),
                "control_standardized_query_threshold_rate": usable.control_query_threshold_rate.mean(),
            }
        )
        layers[layer] = frame[["variant_key", "independent", "query"]]
    both = layers["eqtl"].merge(
        layers["caqtl"], on="variant_key", suffixes=("_eqtl", "_caqtl"), validate="one_to_one"
    )
    both["both_independent"] = both.independent_eqtl & both.independent_caqtl
    save(both, out / "both_independent_locus_background.tsv")
    save(pd.DataFrame(summaries), out / "matched_locus_null.tsv")
    provenance["query_sha256"] = digest(root / "prepared/rfuwas_data1_rfuqtl_grch38.tsv")
    provenance["query_threshold"] = threshold
    (out / "background_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.root, args.out)
