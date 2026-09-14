"""Select chr6/7 Matos summary records for a supplied RFU-QTL variant set.

Offline pilot workflow. Requires pyarrow for the source Parquet tables, outside
scRFU's core dependencies. Verifies released archives before selective extraction.
The archive inventory and all biological outputs belong outside the repository.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import tarfile
from pathlib import Path

import pandas as pd

try:
    from .matos_regulatory_triangulation import adapt_matos_table
except ImportError:
    from matos_regulatory_triangulation import adapt_matos_table

RELEASES = {
    "eqtl": ("18261456", "eQTLs_summary_statistics.tar.gz", "e06e21a30576e6e271d974f922793ea6"),
    "caqtl": ("18317808", "caQTLs_summary_statistics.tar.gz", "0a247f926008e7e7792a7689d327709a"),
}


def digest(path: Path, algorithm: str = "sha256") -> str:
    hasher = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def source_variant_ids(keys: set[str]) -> set[str]:
    """Represent already-normalized GRCh38 keys in the verified Matos ID syntax."""
    result = set()
    for key in keys:
        chrom, pos, ref, alt = key.split(":")
        result.update({f"{chrom}:{pos}[b38]{ref},{alt}", f"chr{chrom}:{pos}[b38]{ref},{alt}"})
    return result


def validate_coverage(filenames: list[str], chromosomes: set[str]) -> None:
    """Missing export coverage is unknown, never a zero-overlap result."""
    for chromosome in sorted(chromosomes):
        local = [name for name in filenames if f".chr{chromosome}." in name]
        if not any(name.endswith(".parquet") for name in local) or not any(
            "independent" in name for name in local
        ):
            raise ValueError(f"Incomplete nominal/conditional source coverage for chr{chromosome}")


def significance_flags(nominal: pd.DataFrame, top: pd.DataFrame | None) -> pd.DataFrame:
    """Apply TensorQTL's phenotype FDR and phenotype-specific nominal threshold.

    Mirrors get_significant_pairs: qval <= .05 and pval_nominal < threshold.
    A row appearing in the misleadingly named sig_cis table is not sufficient.
    """
    if top is None:
        frame = nominal.copy()
        frame["tested"] = True
        frame["nominal_p_lt_005"] = frame.pval_nominal < 0.05
        frame["source_significant"] = pd.array([pd.NA] * len(frame), dtype="boolean")
        frame["source_significance_status"] = "permutation_thresholds_not_released"
        return frame
    columns = ["phenotype_id", "qval", "pval_nominal_threshold", "pval_beta", "pval_perm"]
    if not set(columns).issubset(top):
        raise ValueError("Released permutation table lacks significance metadata")
    if top.phenotype_id.duplicated().any():
        raise ValueError("Ambiguous repeated permutation phenotype")
    metadata = top[columns].rename(columns={c: f"phenotype_{c}" for c in columns[1:]})
    frame = nominal.merge(metadata, on="phenotype_id", how="left", validate="many_to_one")
    if frame.phenotype_qval.isna().any() or frame.phenotype_pval_nominal_threshold.isna().any():
        raise ValueError("Missing phenotype threshold for tested association")
    frame["tested"] = True
    frame["nominal_p_lt_005"] = frame.pval_nominal < 0.05
    frame["source_significant"] = (frame.phenotype_qval <= 0.05) & (
        frame.pval_nominal < frame.phenotype_pval_nominal_threshold
    )
    return frame


def lookup(root: Path, layer: str) -> None:
    import pyarrow.parquet as pq

    record, filename, expected_md5 = RELEASES[layer]
    archive = root / "sources" / "matos" / filename
    prepared = root / "prepared"
    extracted = prepared / f"matos_{layer}_selected"
    checkpoint = prepared / f"matos_{layer}_archive_manifest.json"
    rfu = pd.read_csv(
        prepared / "rfuwas_data1_rfuqtl_grch38.tsv", sep="\t", dtype={"chromosome": str}
    )
    if set(rfu.genome_build) != {"GRCh38"}:
        raise ValueError("This source release requires GRCh38 RFU-QTL inputs")
    chroms = set(rfu.chromosome)
    wanted = source_variant_ids(set(rfu.variant_key))
    if checkpoint.exists():
        manifest = json.loads(checkpoint.read_text())
        if manifest["md5"] != expected_md5 or set(manifest["chromosomes"]) != chroms:
            raise ValueError("Extraction checkpoint differs from requested input")
        for entry in manifest["selected"]:
            if digest(extracted / entry["filename"]) != entry["sha256"]:
                raise ValueError("Selected input differs from verified extraction")
    else:
        print(f"Verifying {filename}", flush=True)
        actual_md5 = digest(archive, "md5")
        if actual_md5 != expected_md5:
            raise ValueError(f"Archive incomplete or checksum mismatch: {actual_md5}")
        manifest = {
            "doi": f"10.5281/zenodo.{record}",
            "archive": filename,
            "bytes": archive.stat().st_size,
            "md5": actual_md5,
            "chromosomes": sorted(chroms),
            "inventory": [],
            "selected": [],
        }
        extracted.mkdir(parents=True, exist_ok=True)
        with tarfile.open(archive, "r|gz") as stream:
            for member in stream:
                manifest["inventory"].append({"name": member.name, "bytes": member.size})
                chrom = re.search(r"\.chr(\d+)\.", member.name)
                if not member.isfile() or chrom is None or chrom[1] not in chroms:
                    continue
                if not (
                    member.name.endswith(".parquet")
                    or "independent" in member.name
                    or "sig_cis" in member.name
                ):
                    continue
                path = extracted / Path(member.name).name
                with stream.extractfile(member) as source, path.open("wb") as out:
                    shutil.copyfileobj(source, out)
                entry = {
                    "filename": path.name,
                    "archive_member": member.name,
                    "bytes": member.size,
                    "sha256": digest(path),
                }
                manifest["selected"].append(entry)
                print(f"Extracted {path.name}", flush=True)
        checkpoint.write_text(json.dumps(manifest, indent=2) + "\n")
    files = [extracted / entry["filename"] for entry in manifest["selected"]]
    validate_coverage([path.name for path in files], chroms)
    nominal, independent, top = [], [], []
    scan_counts = []
    for path in files:
        if path.suffix == ".parquet":
            parquet = pq.ParquetFile(path)
            scan_counts.append({"file": path.name, "rows": parquet.metadata.num_rows})
            retained = []
            for batch in parquet.iter_batches(batch_size=100_000):
                frame = batch.to_pandas()
                frame = frame.loc[frame.variant_id.isin(wanted)].copy()
                if not frame.empty:
                    frame["source_file"] = path.name
                    retained.append(frame)
            nominal.append(pd.concat(retained, ignore_index=True) if retained else frame.iloc[:0])
            print(f"Scanned {path.name}: {sum(map(len, retained))} retained", flush=True)
        else:
            frame = pd.read_csv(path, sep="\t")
            if "phenotype_id" not in frame and "Unnamed: 0" in frame:
                frame = frame.rename(columns={"Unnamed: 0": "phenotype_id"})
            frame = frame.loc[:, ~frame.columns.str.startswith("Unnamed:")]
            frame["source_file"] = path.name
            (independent if "independent" in path.name else top).append(frame)
    if not nominal or not independent:
        raise ValueError("Expected nominal and independent tables; inspect archive inventory")
    raw = significance_flags(
        pd.concat(nominal, ignore_index=True), pd.concat(top, ignore_index=True) if top else None
    )
    indep = pd.concat(independent, ignore_index=True)
    pairs = set(zip(indep.variant_id, indep.phenotype_id, strict=True))
    raw["independent"] = [
        (v, p) in pairs for v, p in zip(raw.variant_id, raw.phenotype_id, strict=True)
    ]
    raw.loc[raw.independent, "source_significant"] = True
    raw.loc[raw.independent, "source_significance_status"] = "released_independent_qtl"
    indep["independent"] = True
    indep["source_significant"] = True
    for name, frame in (
        (f"matos_{layer}_at_rfuqtl_variants.tsv", raw),
        (f"matos_{layer}_independent_chr6_chr7.tsv", indep),
    ):
        adapted = adapt_matos_table(
            frame,
            layer=layer,
            source="Matos CD4 T-cell QTL",
            release=record,
            genome_build="GRCh38",
            format="tensorqtl",
            context="all CD4 T cells",
            allele_order="ref-alt",
        )
        adapted["effect_orientation"] = "unresolved"
        adapted.to_csv(prepared / name, sep="\t", index=False)
    (prepared / f"matos_{layer}_scan_counts.json").write_text(
        json.dumps(scan_counts, indent=2) + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--layer", choices=list(RELEASES), required=True)
    args = parser.parse_args()
    lookup(args.root, args.layer)
