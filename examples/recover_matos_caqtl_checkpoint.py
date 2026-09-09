"""Recover complete members from the known interrupted caQTL archive prefix.

No download and no claim of archive checksum validation. This is a provisional
checkpoint recovery, not a substitute for matos_regulatory_lookup.py on a verified
full archive. A complete full archive should be processed by that script instead.
Run from the repository with its root on PYTHONPATH; requires pyarrow.
"""

import argparse
import json
import shutil
import tarfile
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from examples.matos_regulatory_lookup import digest, source_variant_ids
from examples.matos_regulatory_triangulation import adapt_matos_table


def recover(root: Path) -> None:
    p = root / "prepared"
    a = root / "sources/matos/caQTLs_summary_statistics.tar.gz"
    out = p / "matos_caqtl_partial_selected"
    out.mkdir(exist_ok=True)
    names = {
        "cd4_qsmooth_cpm_chromatin_narrowpeaks.cis_qtl_pairs.chr7.parquet",
        "cd4_qsmooth_cpm_chromatin_narrowpeaks.independent_cis_QTL_pairs.chr6.csv",
    }
    manifest = {
        "archive_complete": False,
        "archive_checksum_verified": False,
        "prefix_bytes": a.stat().st_size,
        "full_archive_expected_bytes": 4273913755,
        "selected": [],
    }
    with tarfile.open(a, "r|gz") as t:
        for m in t:
            name = Path(m.name).name
            if name not in names:
                continue
            path = out / name
            with t.extractfile(m) as source, path.open("wb") as dest:
                shutil.copyfileobj(source, dest)
            assert path.stat().st_size == m.size
            manifest["selected"].append(
                {"archive_member": m.name, "bytes": m.size, "sha256": digest(path)}
            )
            print("Recovered complete member", name, flush=True)
            if len(manifest["selected"]) == len(names):
                break
    r = pd.read_csv(p / "rfuwas_data1_rfuqtl_grch38.tsv", sep="\t")
    wanted = source_variant_ids(set(r.variant_key))
    q = pq.ParquetFile(out / "cd4_qsmooth_cpm_chromatin_narrowpeaks.cis_qtl_pairs.chr7.parquet")
    frames = []
    for b in q.iter_batches(batch_size=100000):
        d = b.to_pandas()
        d = d[d.variant_id.isin(wanted)]
        if len(d):
            frames.append(d)
    d = pd.concat(frames, ignore_index=True)
    d["source_file"] = "cd4_qsmooth_cpm_chromatin_narrowpeaks.cis_qtl_pairs.chr7.parquet"
    d["source_significant"] = pd.NA
    d["independent"] = pd.NA
    for name, df in [
        ("chr7_partial_at_rfuqtl_variants", d),
        (
            "chr6_partial_independent",
            pd.read_csv(
                out / "cd4_qsmooth_cpm_chromatin_narrowpeaks.independent_cis_QTL_pairs.chr6.csv",
                sep="\t",
            ),
        ),
    ]:
        df = adapt_matos_table(
            df,
            layer="caqtl",
            source="Matos CD4 T-cell QTL",
            release="18317808",
            genome_build="GRCh38",
            format="tensorqtl",
            allele_order="ref-alt",
            context="all CD4 T cells",
        )
        df["archive_checksum_verified"] = False
        df.to_csv(p / f"matos_caqtl_{name}.tsv", sep="\t", index=False)
    manifest["chr7_nominal_rows_scanned"] = q.metadata.num_rows
    manifest["chr7_nominal_rows_retained"] = len(d)
    (p / "matos_caqtl_partial_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    recover(parser.parse_args().root)
