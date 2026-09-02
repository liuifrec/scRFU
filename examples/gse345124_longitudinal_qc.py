#!/usr/bin/env python3
"""Audit processed GSE345124 TCR tables without running RFU assignment."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

import scrfu

NAME_PATTERN = re.compile(
    r"(?P<accession>GSM\d+)_(?P<donor>DCP\d+)_(?P<time>D\d+)_(?P<compartment>CD[48])_PBMC_TCR.tsv.gz"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_state() -> tuple[str | None, bool | None]:
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True
        ).stdout.strip()
        dirty = bool(
            subprocess.run(
                ["git", "status", "--porcelain"],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
        )
        return commit, dirty
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None, None


def run(args: argparse.Namespace) -> dict[str, Any]:
    source = args.input_dir.expanduser().resolve()
    archive = args.archive.expanduser().resolve()
    output = args.output_dir.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    files = sorted(source.glob("*.tsv.gz"))
    if not files:
        raise FileNotFoundError(f"No processed TSV.GZ files found in {source}")

    global_eligible: set[str] = set()
    sample_rows: list[dict[str, Any]] = []
    schema_counts: dict[str, int] = {}
    for path in files:
        match = NAME_PATTERN.fullmatch(path.name)
        if match is None:
            raise ValueError(f"Unexpected GSE345124 filename: {path.name}")
        columns = pd.read_csv(path, sep="\t", nrows=0).columns
        legacy = "aminoAcid" in columns
        schema = "adaptive_legacy" if legacy else "adaptive_current"
        schema_counts[schema] = schema_counts.get(schema, 0) + 1
        aa_col = "aminoAcid" if legacy else "amino_acid"
        v_col = "vGeneName" if legacy else "v_gene"
        status_col = "sequenceStatus" if legacy else "frame_type"
        count_col = "count (templates/reads)" if legacy else "templates"
        frame = pd.read_csv(path, sep="\t", usecols=[aa_col, v_col, status_col, count_col])
        productive = frame[status_col].eq("In")
        amino_acids = frame.loc[productive, aa_col].astype("string")
        eligible = amino_acids.str.startswith("C", na=False) & amino_acids.str.endswith(
            ("F", "W"), na=False
        )
        eligible_values = amino_acids.loc[eligible].dropna()
        global_eligible.update(eligible_values.tolist())
        identifiers = match.groupdict()
        sample_rows.append(
            {
                **identifiers,
                "sample_id": path.name.removesuffix("_TCR.tsv.gz"),
                "schema": schema,
                "input_rows": len(frame),
                "in_frame_rows": int(productive.sum()),
                "eligible_rows": int(eligible.sum()),
                "unique_eligible_cdr3": int(eligible_values.nunique()),
                "missing_v_in_frame": int(frame.loc[productive, v_col].isna().sum()),
                "template_count_in_frame": float(
                    pd.to_numeric(frame.loc[productive, count_col], errors="coerce").sum()
                ),
                "file_size_bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
        )
    sample_qc = pd.DataFrame(sample_rows)
    if sample_qc["sample_id"].duplicated().any():
        raise RuntimeError("GSE345124 sample identifiers are not unique.")
    donor_visits = sample_qc[["donor", "time"]].drop_duplicates()
    visits_per_donor = donor_visits.groupby("donor")["time"].nunique()
    expected_compartments = sample_qc.groupby(["donor", "time"])["compartment"].nunique()
    git_commit, git_dirty = _git_state()
    sample_path = output / "sample_qc.tsv"
    sample_qc.to_csv(sample_path, sep="\t", index=False)
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "GSE345124 processed bulk TCR longitudinal candidate",
        "public_source": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE345124",
        "download_source": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE345nnn/GSE345124/suppl/GSE345124_RAW.tar",
        "archive": {
            "filename": archive.name,
            "size_bytes": archive.stat().st_size,
            "sha256": _sha256(archive),
        },
        "software": {
            "scrfu_version": scrfu.__version__,
            "python_version": platform.python_version(),
            "git_commit": git_commit,
            "git_dirty": git_dirty,
        },
        "counts": {
            "files": len(sample_qc),
            "schemas": schema_counts,
            "donors": int(sample_qc["donor"].nunique()),
            "donor_timepoints": len(donor_visits),
            "donors_with_three_timepoints": int(visits_per_donor.eq(3).sum()),
            "donors_with_two_timepoints": int(visits_per_donor.eq(2).sum()),
            "timepoint_file_counts": {
                str(key): int(value) for key, value in sample_qc["time"].value_counts().items()
            },
            "compartment_file_counts": {
                str(key): int(value)
                for key, value in sample_qc["compartment"].value_counts().items()
            },
            "donor_timepoints_with_both_compartments": int(expected_compartments.eq(2).sum()),
            "input_rows": int(sample_qc["input_rows"].sum()),
            "in_frame_rows": int(sample_qc["in_frame_rows"].sum()),
            "eligible_rows": int(sample_qc["eligible_rows"].sum()),
            "unique_eligible_cdr3": len(global_eligible),
            "missing_v_in_frame": int(sample_qc["missing_v_in_frame"].sum()),
        },
        "execution_gate": {
            "status": "deferred",
            "reason": (
                "A full frozen-reference run would require 2,258,994 distinct eligible CDR3 "
                "queries and genuine Scirpy comparison over 3,038,924 in-frame rows; this is "
                "outside the bounded optional-validation gate and was not silently subsampled."
            ),
        },
        "outputs": {
            sample_path.name: {
                "rows": len(sample_qc),
                "size_bytes": sample_path.stat().st_size,
                "sha256": _sha256(sample_path),
            }
        },
    }
    (output / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing extracted public GSE345124 *_TCR.tsv.gz files.",
    )
    parser.add_argument(
        "--archive",
        type=Path,
        required=True,
        help="Downloaded public GSE345124_RAW.tar archive to hash for provenance.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Explicit external directory for QC tables and the run manifest.",
    )
    args = parser.parse_args()
    print(json.dumps(run(args), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
