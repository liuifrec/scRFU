"""Audit and reuse the restored RP1-14 assignments; biological outputs stay external.

The historical per-receptor export truncated input labels after EncodeRepertoire
filtered non-C-starting sequences. This loader repairs *alignment*, never RFU
values, only after an exact source-prefix/length reconciliation. It fails closed
on any other layout. No mixed-cohort matrix or newly assigned full repertoire is
accepted. A small, explicitly saved parity assay is the only RFU recomputation.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import inspect
import itertools
import json
import subprocess
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

import scrfu
from scrfu.io import file_sha256

from .radiation_methods import REPO, completed, json_save, save, validate_boundary

DATASET = "RP1-14_published_authorized"
CODONS = dict(
    zip(
        map("".join, itertools.product("TCAG", repeat=3)),
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        strict=True,
    )
)
RAW_COLUMNS = {
    "nucleotide",
    "aminoAcid",
    "count",
    "count (reads)",
    "frequencyCount (%)",
    "cdr3Length",
    "vIndex",
    "vMaxResolved",
    "jMaxResolved",
    "vGeneName",
    "jGeneName",
    "sequenceStatus",
}


def registry_from_metadata(metadata: pd.DataFrame) -> pd.DataFrame:
    """Order visits by measured collection age, never by file suffix."""
    required = ["Sample ID", "CD4 or CD8", "Age at collection"]
    if not set(required).issubset(metadata):
        raise ValueError("Source sample/compartment/collection-age metadata required.")
    frame = metadata.dropna(subset=["Sample ID"])[required].copy()
    frame.columns = ["source_sample", "compartment", "age_years"]
    parsed = frame.source_sample.str.extract(
        r"^D(?P<source_donor>\d+)-(?P<parsed_compartment>CD[48])_(?P<source_visit>\d+)$"
    )
    if parsed.isna().any().any() or not parsed.parsed_compartment.eq(frame.compartment).all():
        raise ValueError("Invalid published sample naming or conflicting compartment.")
    frame = frame.join(parsed.drop(columns="parsed_compartment"))
    if frame.source_sample.duplicated().any():
        raise ValueError("Duplicated sample metadata.")
    frame["age_years"] = pd.to_numeric(frame.age_years, errors="raise")
    if not np.isfinite(frame.age_years).all():
        raise ValueError("Missing or nonfinite chronological metadata.")
    if frame.groupby(["source_donor", "source_visit"]).age_years.nunique().gt(1).any():
        raise ValueError("CD4/CD8 collection ages disagree within a source visit.")
    visits = frame[["source_donor", "source_visit", "age_years"]].drop_duplicates()
    if visits.duplicated(["source_donor", "age_years"]).any():
        raise ValueError("Multiple visit suffixes have the same collection age.")
    donors = sorted(frame.source_donor.unique(), key=int)
    frame["donor"] = frame.source_donor.map({x: f"RP{i + 1:02}" for i, x in enumerate(donors)})
    frame["visit"] = frame.groupby("donor").age_years.rank(method="dense").astype(int)
    frame["elapsed_years"] = frame.age_years - frame.groupby("donor").age_years.transform("min")
    frame["sample_id"] = frame.donor + "_" + frame.compartment + "_V" + frame.visit.astype(str)
    frame["dataset_id"], frame["authorized"] = DATASET, True
    validate_boundary(frame, frame)
    return frame.sort_values(["donor", "compartment", "visit"]).reset_index(drop=True)


def normalize_raw(frame: pd.DataFrame) -> pd.DataFrame:
    """Accept either documented read-count header; never infer counts from percent."""
    found = [x for x in ("count", "count (reads)") if x in frame]
    if len(found) != 1:
        raise ValueError(
            "Exactly one source read-count column is required; transformed input rejected."
        )
    result = frame.rename(columns={found[0]: "read_count"}).copy()
    counts = pd.to_numeric(result.read_count, errors="raise")
    if not np.isfinite(counts).all() or counts.lt(0).any() or not counts.eq(np.floor(counts)).all():
        raise ValueError("Source reads must be nonnegative integer counts, not transformed mass.")
    if counts.sum() <= 0:
        raise ValueError("Empty source read denominator.")
    result["read_count"] = counts.astype("int64")
    if "frequencyCount (%)" in result:
        expected = counts / counts.sum() * 100
        if not np.allclose(result["frequencyCount (%)"], expected, atol=5e-6, rtol=0):
            raise ValueError("Source percentage and full read denominator disagree.")
    return result


def restore_alignment(
    raw: pd.DataFrame, saved: pd.DataFrame, top_n: int
) -> tuple[pd.DataFrame, dict]:
    """Recover input-row identity without altering saved labels, scores or flags."""
    if not raw.read_count.is_monotonic_decreasing:
        raise ValueError("Historical source order is not decreasing read count.")
    selected = raw[raw.aminoAcid.notna()].head(top_n).copy()
    selected["source_row"] = selected.index
    selected = selected.reset_index(drop=True)
    encoded = selected[selected.aminoAcid.str.startswith("C")].reset_index(drop=True)
    if len(saved) != len(encoded):
        raise ValueError("Saved assignment count does not match upstream C-start filter.")
    prefix = selected.iloc[: len(saved)].reset_index(drop=True)
    if not saved.cdr3_aa.eq(prefix.aminoAcid).all() or not saved.trbv.eq(prefix.vMaxResolved).all():
        raise ValueError("Historical prefix alignment could not be verified; do not guess a join.")
    labels = pd.to_numeric(saved.rfu, errors="raise")
    nonmissing = labels.dropna()
    if not nonmissing.between(1, 5000).all() or not nonmissing.eq(np.floor(nonmissing)).all():
        raise ValueError("Historical one-based RFU labels are invalid.")
    scores = pd.to_numeric(saved.max_cor, errors="raise")
    flags = saved.pass_thr.astype("string").str.upper().map({"TRUE": True, "FALSE": False})
    if not flags[flags.notna()].eq(scores[flags.notna()].ge(0.6)).all():
        raise ValueError("Historical qualification disagrees with threshold 0.6.")
    encoded["rfu_numeric"] = labels.astype("Int64")
    encoded["rfu"] = "RFU" + labels.astype("Int64").astype("string")
    encoded["max_cor"] = scores
    encoded["pass_threshold"] = flags.fillna(False).astype(bool)
    return encoded, {
        "selected_translated_rows": len(selected),
        "saved_rows": len(saved),
        "non_C_rows_excluded_by_upstream": len(selected) - len(encoded),
        "mislabelled_aa_rows_repaired": int(saved.cdr3_aa.ne(encoded.aminoAcid).sum()),
        "missing_rfu": int(labels.isna().sum()),
        "missing_saved_frequency": int(saved.freq.isna().sum()),
        "repair": "assignment_vector_unchanged_reattached_to_upstream_filtered_source_order",
    }


def receptor_features(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    nt = [
        s[int(i) : int(i) + int(n)]
        for s, i, n in zip(result.nucleotide, result.vIndex, result.cdr3Length, strict=True)
    ]
    translated = ["".join(CODONS.get(s[i : i + 3], "X") for i in range(0, len(s), 3)) for s in nt]
    if not all(len(s) == n for s, n in zip(nt, result.cdr3Length, strict=True)):
        raise ValueError("CDR3 nucleotide interval extends beyond the observed sequence.")
    if not np.array_equal(np.asarray(translated), result.aminoAcid.to_numpy()):
        raise ValueError("Nucleotide CDR3 translation does not match source amino acids.")
    result["cdr3nt"] = nt
    result["clone"] = (
        result.cdr3nt
        + "|"
        + result.vMaxResolved.fillna("unresolved")
        + "|"
        + result.jMaxResolved.fillna("unresolved")
    )
    for segment in ("v", "j"):
        gene = result[f"{segment}GeneName"].fillna("unresolved").astype(str)
        fallback = (
            result[f"{segment}MaxResolved"]
            .fillna("unresolved")
            .str.replace(r"\*.*$", "", regex=True)
        )
        resolved = ~gene.str.lower().eq("unresolved")
        result[f"{segment}_resolved"] = resolved
        result[f"{segment}_call"] = gene.where(resolved, "unresolved_or_family:" + fallback)
    result["vj"] = result.v_call + "|" + result.j_call
    result["length_bin"] = result.aminoAcid.str.len() // 5
    result["productive_canonical"] = result.sequenceStatus.eq(
        "In"
    ) & result.aminoAcid.str.fullmatch(r"C[ACDEFGHIKLMNPQRSTVWY]+")
    return result


def assert_fixed_mapping(frame: pd.DataFrame) -> None:
    for keys, values in [
        ("aminoAcid", ["rfu", "max_cor"]),
        ("clone", ["rfu", "v_call", "j_call", "length_bin"]),
    ]:
        if frame.groupby(keys)[values].nunique().gt(1).any().any():
            raise ValueError("Reused receptors do not have a fixed assignment/feature mapping.")


def reconcile_matrix(matrix: pd.DataFrame, assigned: dict[str, pd.DataFrame], scale: int) -> float:
    if len(matrix) != 5000 or set(matrix.columns) != set(assigned):
        raise ValueError(
            "Historical matrix must contain exactly the authorized saved samples and 5000 rows."
        )
    max_error = 0.0
    for sample, frame in assigned.items():
        expected = (
            frame.rfu_numeric.value_counts().reindex(range(1, 5001), fill_value=0)
            / len(frame)
            * scale
        )
        max_error = max(
            max_error, float(np.max(np.abs(expected.to_numpy() - matrix[sample].to_numpy())))
        )
    if max_error > 1e-8:
        raise ValueError(
            "Matrix is transformed, weighted differently or not the saved assignment matrix."
        )
    return max_error


def parity_check(frame: pd.DataFrame, reference: Path, out: Path, config: dict) -> dict:
    """Bounded compatibility check; never used to replace historical assignments."""
    eligible = frame[frame.productive_canonical & frame.rfu.notna()].copy()
    parts = []
    for i, (_, group) in enumerate(eligible.groupby("sample_id", sort=True)):
        parts.append(
            group.drop_duplicates("aminoAcid").sample(
                min(config["parity_receptors_per_sample"], group.aminoAcid.nunique()),
                random_state=config["random_seed"] + i,
            )
        )
    sample = pd.concat(parts).drop_duplicates("aminoAcid").reset_index(drop=True)
    sample["cell_id"] = [f"parity_{i:04}" for i in range(len(sample))]
    receptors = pd.DataFrame(
        {
            "cell_id": sample.cell_id,
            "chain": "TRB",
            "cdr3aa": sample.aminoAcid,
            "v_call": sample.v_call,
            "source_adapter": "rp1_14_bounded_parity",
        }
    )
    cache = out / "bounded_parity"
    pin_path = cache / "parity_integrity.json"
    expected = {
        "reference": {
            name: file_sha256(reference / name)
            for name in ["km5000noMax.Rdata", "RFU.R", "trimerMDSfit_small.Rdata"]
        },
        "wrapper_sha256": file_sha256(REPO / "r/run_rfu_repo.R"),
        "threshold": config["threshold"],
    }
    if pin_path.exists():
        pin = json.loads(pin_path.read_text())
        if any(pin[k] != v for k, v in expected.items()) or any(
            file_sha256(cache / name) != value for name, value in pin["files"].items()
        ):
            raise ValueError("Bounded parity cache scientific inputs or outputs changed.")
        inputs = pd.read_csv(cache / "rfu_in.tsv", sep="\t")
        outputs = pd.read_csv(cache / "rfu_out.tsv", sep="\t")
        if set(inputs.cdr3aa) != set(receptors.cdr3aa) or not inputs.set_index(
            "cdr3aa"
        ).trbv.sort_index().equals(
            receptors.set_index("cdr3aa").v_call.sort_index().rename("trbv")
        ):
            raise ValueError("Bounded parity receptor set changed.")
        values = inputs.merge(outputs, on="unique_sequence_id", validate="one_to_one")
        values = receptors[["cell_id", "cdr3aa"]].merge(values, on="cdr3aa", validate="one_to_one")
        values = values.rename(columns={"rfu_label": "rfu_label_nearest"})
    else:
        result = scrfu.tl.assign_rfu(
            receptors,
            rfu_dir=reference,
            mode="standard",
            threshold=config["threshold"],
            rscript_bin="/usr/bin/Rscript",
            workdir=cache,
            max_workers=1,
            wrapper_r_path=REPO / "r/run_rfu_repo.R",
        )
        values = result.per_row
        json_save(
            {
                **expected,
                "files": {
                    name: file_sha256(cache / name) for name in ["rfu_in.tsv", "rfu_out.tsv"]
                },
            },
            pin_path,
        )
    compare = sample[["cell_id", "sample_id", "aminoAcid", "rfu", "max_cor"]].merge(
        values[["cell_id", "rfu_label_nearest", "rfu_score"]].rename(
            columns={"rfu_score": "current_max_cor"}
        ),
        on="cell_id",
        validate="one_to_one",
    )
    compare["same_label"] = compare.rfu.eq(compare.rfu_label_nearest)
    compare["score_error"] = (compare.max_cor - compare.current_max_cor).abs()
    save(compare, out, "bounded_reference_parity.tsv")
    return {
        "n": len(compare),
        "matching_labels": int(compare.same_label.sum()),
        "max_score_error": float(compare.score_error.max()),
        "interpretation": "bounded computational parity; historical execution/reference hash not independently recorded",
    }


def run(args: argparse.Namespace) -> None:
    config = json.loads(args.config.read_text())
    if config["dataset_id"] != DATASET:
        raise ValueError("Only the published RP1-14 authorization scope is accepted.")
    out = args.out.resolve()
    if out.is_relative_to(REPO):
        raise ValueError("Biological prepared outputs must be outside the repository.")
    inputs = sorted([*args.raw.glob("*.xlsx"), *args.raw.glob("*.zip"), *args.rfu.rglob("*.tsv")])
    hashes = {str(p.resolve()): file_sha256(p) for p in inputs}
    identity = {
        "inputs": hashes,
        "config": file_sha256(args.config),
        "script": file_sha256(Path(__file__)),
        "reference": {
            p: file_sha256(args.reference / p)
            for p in ("km5000noMax.Rdata", "RFU.R", "trimerMDSfit_small.Rdata")
        },
    }
    if identity["reference"]["km5000noMax.Rdata"] != config["reference_sha256"]:
        raise ValueError("Frozen reference hash mismatch.")
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    if completed(out, fingerprint):
        print("RP1-14 preparation verified and reused without recomputation.")
        return
    out.mkdir(parents=True, exist_ok=True)
    (out / "full_productive_counts").mkdir(exist_ok=True)
    json_save({"status": "running", "fingerprint": fingerprint}, out / "completion.json")
    preparation_identity = {
        "inputs": hashes,
        "config": identity["config"],
        "functions": {
            name: hashlib.sha256(inspect.getsource(globals()[name]).encode()).hexdigest()
            for name in [
                "registry_from_metadata",
                "normalize_raw",
                "restore_alignment",
                "receptor_features",
                "assert_fixed_mapping",
                "reconcile_matrix",
            ]
        },
    }
    stage_path = out / "prepared_stage.json"
    if stage_path.exists():
        stage = json.loads(stage_path.read_text())
        if stage.get("status") != "complete" or stage["identity"] != preparation_identity:
            raise ValueError("Prepared stage scientific inputs changed; use a new directory.")
        for name, checksum in stage["outputs"].items():
            if file_sha256(out / name) != checksum:
                raise ValueError("Prepared stage output checksum failed.")
        frame = pd.read_parquet(out / "reused_receptor_assignments.parquet")
        registry = pd.read_csv(out / "source_sample_crosswalk.tsv", sep="\t")
        coverage = pd.read_csv(out / "sample_registry.tsv", sep="\t").to_dict(orient="records")
        audits = pd.read_csv(out / "historical_alignment_repair.tsv", sep="\t").to_dict(
            orient="records"
        )
        error = stage["matrix_error"]
        print("Verified preparation stage reused", flush=True)
    else:
        registry = registry_from_metadata(pd.read_excel(args.raw / "TCR-seq gender and age.xlsx"))
        lookup = {p.name.replace(".cdr3_rfu.tsv", ""): p for p in args.rfu.rglob("*.cdr3_rfu.tsv")}
        if not set(lookup).issubset(registry.source_sample):
            raise ValueError(
                "Assignment contains a sample outside the authorized published registry."
            )
        audits, coverage, frames, assigned, inventory = [], [], [], {}, []
        for path in inputs:
            inventory.append(
                {
                    "current_path": str(path.resolve()),
                    "dataset": DATASET,
                    "bytes": path.stat().st_size,
                    "sha256": hashes[str(path.resolve())],
                    "format": path.suffix,
                    "status": "authorized_analysis_not_for_git_redistribution",
                }
            )
        seen = set()
        for archive in sorted(args.raw.glob("*.zip")):
            with zipfile.ZipFile(archive) as zf:
                for member in zf.namelist():
                    if not member.endswith(".tsv"):
                        continue
                    source_sample = Path(member).stem
                    if source_sample in seen or source_sample not in set(registry.source_sample):
                        raise ValueError("Repeated or unauthorized source sample.")
                    seen.add(source_sample)
                    meta = registry.set_index("source_sample").loc[source_sample].to_dict()
                    raw = normalize_raw(
                        pd.read_csv(zf.open(member), sep="\t", usecols=lambda x: x in RAW_COLUMNS)
                    )
                    full = receptor_features(
                        raw[raw.sequenceStatus.eq("In") & raw.aminoAcid.notna()]
                    )
                    full_counts = full.groupby(
                        ["clone", "v_call", "j_call"], as_index=False
                    ).read_count.sum()
                    full_counts.to_parquet(
                        out / "full_productive_counts" / f"{meta['sample_id']}.parquet", index=False
                    )
                    if source_sample not in lookup:
                        coverage.append(
                            {
                                **meta,
                                "status": "missing_saved_assignment",
                                "source_read_depth": int(raw.read_count.sum()),
                            }
                        )
                        continue
                    saved = pd.read_csv(lookup[source_sample], sep="\t")
                    aligned, audit = restore_alignment(
                        raw, saved, config["historical_top_translated_rows"]
                    )
                    aligned = receptor_features(aligned)
                    assigned[source_sample] = aligned
                    aligned = aligned.assign(**meta)
                    frames.append(aligned)
                    audits.append({"sample_id": meta["sample_id"], **audit})
                    primary = aligned[
                        aligned.productive_canonical & aligned.pass_threshold & aligned.rfu.notna()
                    ]
                    coverage.append(
                        {
                            **meta,
                            "status": "analyzed",
                            "source_rows": len(raw),
                            "source_read_depth": int(raw.read_count.sum()),
                            "productive_rows": len(full),
                            "productive_read_depth": int(full.read_count.sum()),
                            "productive_unique_receptors": full.clone.nunique(),
                            "historic_assigned_rows": len(aligned),
                            "historic_assigned_reads": int(aligned.read_count.sum()),
                            "primary_rows": len(primary),
                            "primary_unique_receptors": primary.clone.nunique(),
                            "primary_reads": int(primary.read_count.sum()),
                            "observed_primary_rfus": primary.rfu.nunique(),
                            "primary_fraction_all_reads": float(
                                primary.read_count.sum() / raw.read_count.sum()
                            ),
                            "primary_fraction_productive_reads": float(
                                primary.read_count.sum() / full.read_count.sum()
                            ),
                            "resolved_V_fraction_primary_reads": float(
                                primary.loc[primary.v_resolved, "read_count"].sum()
                                / primary.read_count.sum()
                            ),
                            "resolved_J_fraction_primary_reads": float(
                                primary.loc[primary.j_resolved, "read_count"].sum()
                                / primary.read_count.sum()
                            ),
                        }
                    )
                    print(
                        f"Prepared {meta['sample_id']}: {len(raw)} source rows; {len(primary)} qualified historical rows",
                        flush=True,
                    )
        frame = pd.concat(frames, ignore_index=True)
        assert_fixed_mapping(frame)
        error = reconcile_matrix(
            pd.read_csv(args.rfu / "RFU_matrix.tsv", sep="\t"),
            assigned,
            config["historical_matrix_scale"],
        )
        for row in registry[~registry.source_sample.isin(seen)].to_dict(orient="records"):
            coverage.append({**row, "status": "missing_raw_sample"})
        frame.to_parquet(out / "reused_receptor_assignments.parquet", index=False)
        save(pd.DataFrame(coverage), out, "sample_registry.tsv")
        save(registry, out, "source_sample_crosswalk.tsv")
        save(pd.DataFrame(audits), out, "historical_alignment_repair.tsv")
        save(pd.DataFrame(inventory), out, "rp1_14_asset_manifest.tsv")
        names = [
            "reused_receptor_assignments.parquet",
            "source_sample_crosswalk.tsv",
            "sample_registry.tsv",
            "historical_alignment_repair.tsv",
            "rp1_14_asset_manifest.tsv",
        ]
        paths = [out / n for n in names] + list((out / "full_productive_counts").glob("*.parquet"))
        json_save(
            {
                "status": "complete",
                "identity": preparation_identity,
                "matrix_error": error,
                "outputs": {str(p.relative_to(out)): file_sha256(p) for p in paths},
            },
            stage_path,
        )
    parity = parity_check(frame, args.reference, out, config)
    if parity["matching_labels"] != parity["n"] or parity["max_score_error"] > 1e-8:
        raise ValueError("Bounded current-reference parity failed; inspect before analysis.")
    summary = {
        "samples": len(coverage),
        "donors": registry.donor.nunique(),
        "historic_rows": len(frame),
        "historical_matrix_max_error": error,
        "mislabelled_aa_rows_repaired": int(
            pd.DataFrame(audits).mislabelled_aa_rows_repaired.sum()
        ),
        "repaired_assignment_conflicts": 0,
        "bounded_parity": parity,
        "cross_study_label_identity": "consistent_with_bounded_parity_not_full_historical_execution_provenance",
        "matrix_semantics": "unweighted saved assignment rows per 10000; nearest labels including below threshold; not read counts",
        "source_abundance_unit": "sequencing_reads_not_cells_or_observed_templates",
        "source_reference_provenance": "historical reference/checksum and generating wrapper unavailable",
    }
    json_save(summary, out / "evidence_counts.json")
    json_save(
        {
            **identity,
            "fingerprint": fingerprint,
            "command": sys.argv,
            "python": sys.executable,
            "scrfu_import": scrfu.__file__,
            "git_sha": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
            ).strip(),
            "rscript": "/usr/bin/Rscript",
            "r_version": subprocess.check_output(
                ["/usr/bin/Rscript", "--version"], text=True
            ).strip(),
            "versions": {
                k: importlib.metadata.version(k)
                for k in ("pandas", "numpy", "pyarrow", "openpyxl", "scrfu")
            },
            "summary": summary,
        },
        out / "rp1_14_provenance.json",
    )
    outputs = {
        str(p.relative_to(out)): file_sha256(p)
        for p in sorted(out.rglob("*"))
        if p.is_file() and p.name != "completion.json"
    }
    json_save(
        {"status": "complete", "fingerprint": fingerprint, "outputs": outputs},
        out / "completion.json",
    )
    print(json.dumps(summary, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--rfu", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--reference", type=Path, default=Path("/home/liuyuchen/ext/RFU-official"))
    parser.add_argument("--config", type=Path, default=REPO / "manuscript/config/rp1_14_v1.json")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
