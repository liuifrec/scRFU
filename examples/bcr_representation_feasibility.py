"""Evaluate exact, outcome-independent BCR representation candidates without clustering."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from scrfu import __version__


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _demux_map(directory: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for path in sorted(directory.glob("*.best")):
        sample = path.name.removesuffix("_demux.best")
        frame = pd.read_csv(path, sep="\t", low_memory=False)
        single = frame["DROPLET.TYPE"].astype(str).eq("SNG")
        for barcode, donor in frame.loc[single, ["BARCODE", "SNG.BEST.GUESS"]].itertuples(
            index=False, name=None
        ):
            result[f"{sample}::{barcode}"] = str(donor)
    return result


def _normalize_cell(value: str) -> str:
    return re.sub(r"^GSM\d+_", "", value)


def _key(frame: pd.DataFrame, fields: list[str]) -> pd.Series:
    available = frame[fields].notna().all(axis=1)
    result = pd.Series(pd.NA, index=frame.index, dtype="string")
    result.loc[available] = frame.loc[available, fields].astype(str).agg("|".join, axis=1)
    return result


def _weighted_group_purity(frame: pd.DataFrame, key: str, label: str) -> float | None:
    working = frame.loc[frame[key].notna() & frame[label].notna(), [key, label]]
    if working.empty:
        return None
    counts = working.groupby([key, label], observed=True).size().rename("n").reset_index()
    totals = counts.groupby(key, observed=True)["n"].sum()
    maxima = counts.groupby(key, observed=True)["n"].max()
    return float(maxima.sum() / totals.sum())


def _candidate(
    development: pd.DataFrame,
    validation: pd.DataFrame,
    *,
    label: str,
    fields: list[str],
    donor_col: str,
    seeds: list[int],
) -> dict[str, Any]:
    dev = development.copy()
    val = validation.copy()
    dev["candidate_key"] = _key(dev, fields)
    val["candidate_key"] = _key(val, fields)
    eligible_dev = dev.loc[dev["candidate_key"].notna()].copy()
    eligible_val = val.loc[val["candidate_key"].notna()].copy()
    groups = eligible_dev["candidate_key"].value_counts()
    reference = set(groups.index.astype(str))
    mapped = eligible_val["candidate_key"].astype(str).isin(reference)
    downsampling: list[float] = []
    for seed in seeds:
        selected = eligible_dev.sample(frac=0.5, random_state=seed)
        heldout = eligible_dev.drop(selected.index)
        downsampling.append(
            float(heldout["candidate_key"].isin(set(selected["candidate_key"])).mean())
            if len(heldout)
            else float("nan")
        )
    donor_purity = _weighted_group_purity(eligible_dev, "candidate_key", donor_col)
    return {
        "candidate": label,
        "fields": fields,
        "distance": "exact equality on all available required fields",
        "missing_data_behavior": "unassigned when any required field is missing",
        "development_eligible_cells": len(eligible_dev),
        "development_eligible_fraction": len(eligible_dev) / max(1, len(dev)),
        "development_group_count": len(groups),
        "development_singleton_group_fraction": float(groups.eq(1).mean()) if len(groups) else None,
        "development_largest_group_fraction": float(groups.max() / len(eligible_dev))
        if len(groups)
        else None,
        "development_donor_assignment_fraction": float(eligible_dev[donor_col].notna().mean()),
        "development_weighted_donor_purity": donor_purity,
        "development_weighted_clonotype_purity": _weighted_group_purity(
            eligible_dev, "candidate_key", "heavy_clonotype_id"
        ),
        "development_weighted_isotype_purity": _weighted_group_purity(
            eligible_dev, "candidate_key", "heavy_isotype"
        ),
        "half_sample_holdout_coverage_mean": float(np.nanmean(downsampling)),
        "half_sample_holdout_coverage_min": float(np.nanmin(downsampling)),
        "validation_eligible_cells": len(eligible_val),
        "validation_eligible_fraction": len(eligible_val) / max(1, len(val)),
        "frozen_exact_reference_coverage": float(mapped.mean()) if len(mapped) else None,
        "frozen_exact_reference_unassigned_fraction": float((~mapped).mean())
        if len(mapped)
        else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--development-features", type=Path, required=True)
    parser.add_argument("--validation-features", type=Path, required=True)
    parser.add_argument("--development-demux", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    development = pd.read_csv(args.development_features, sep="\t", low_memory=False)
    validation = pd.read_csv(args.validation_features, sep="\t", low_memory=False)
    donor_map = _demux_map(args.development_demux)
    development["_donor"] = (
        development["cell_id"].astype(str).map(lambda value: donor_map.get(_normalize_cell(value)))
    )
    seeds = [20260825, 20260826, 20260827, 20260828, 20260829]
    candidates = [
        _candidate(
            development,
            validation,
            label="A_heavy_cdr3_exact",
            fields=["heavy_cdr3aa"],
            donor_col="_donor",
            seeds=seeds,
        ),
        _candidate(
            development,
            validation,
            label="B_heavy_vj_cdr3_exact",
            fields=["heavy_v_call", "heavy_j_call", "heavy_cdr3aa"],
            donor_col="_donor",
            seeds=seeds,
        ),
        _candidate(
            development,
            validation,
            label="C_paired_heavy_light_exact",
            fields=[
                "heavy_v_call",
                "heavy_j_call",
                "heavy_cdr3aa",
                "light_v_call",
                "light_j_call",
                "light_cdr3aa",
            ],
            donor_col="_donor",
            seeds=seeds,
        ),
    ]
    report = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "scrfu_version": __version__,
        "development_dataset": "GSE219098",
        "validation_dataset": "GSE266519",
        "development_feature_sha256": _sha256(args.development_features),
        "validation_feature_sha256": _sha256(args.validation_features),
        "random_seeds": seeds,
        "outcome_labels_used_for_construction": False,
        "candidates": candidates,
        "unavailable_candidates": [
            {
                "candidate": "D_maturation_aware",
                "reason": "mutation frequency and germline identity are absent in both acquired receptor tables",
            },
            {
                "candidate": "E_clonal_family_aware",
                "reason": "inferred clonal-family identifiers are absent in both acquired receptor tables",
            },
        ],
        "decision": "NO-GO",
        "major_failures": [
            "No candidate supplies a justified non-exact receptor distance without arbitrary weights.",
            "The acquired tables do not expose SHM/germline identity or inferred clonal families.",
            "Exact candidates are baselines, not novel receptor-state groups, and have poor frozen cross-dataset coverage.",
            "The development cohort has only three donors and incomplete donor assignment; the validation aggregate lacks compatible donor metadata.",
        ],
        "functional_reference_constructed": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
