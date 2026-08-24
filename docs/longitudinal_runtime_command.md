# Governed longitudinal runtime command

Checked on 2026-08-24. `SCRFU_LONGITUDINAL_RECEPTORS`,
`SCRFU_LONGITUDINAL_METADATA`, and `SCRFU_LONGITUDINAL_OUTDIR` are unset. No
private-file search was performed and no longitudinal result was simulated.

The input contract is deliberately explicit:

- receptors: CSV/TSV, optionally gzip-compressed, satisfying the canonical
  receptor schema and containing the declared biological-sample key;
- metadata: CSV/TSV, optionally gzip-compressed, with one nonconflicting mapping
  from sample to donor, numeric or ordered time, and optional compartment;
- sample, donor, and time key names: required environment variables for the
  governed run; no silent default;
- compartment key: optional, but required for compartment analyses;
- output: an external governed directory, never the repository.

## Exact preflight command

Set the paths and key names, then run this command from the repository. It writes
only a de-identified count/hash preflight report to the external output
directory. It does not print participant or sample identifiers.

```bash
export SCRFU_LONGITUDINAL_RECEPTORS=/path/to/private/receptors.tsv.gz
export SCRFU_LONGITUDINAL_METADATA=/path/to/private/metadata.tsv.gz
export SCRFU_LONGITUDINAL_OUTDIR=/path/to/private/scRFU-output
export SCRFU_LONGITUDINAL_SAMPLE_KEY=sample_id
export SCRFU_LONGITUDINAL_DONOR_KEY=donor_id
export SCRFU_LONGITUDINAL_TIME_KEY=time
export SCRFU_LONGITUDINAL_COMPARTMENT_KEY=compartment

python - <<'PY'
import hashlib
import json
import os
from pathlib import Path

import pandas as pd
import scrfu

required = [
    "SCRFU_LONGITUDINAL_RECEPTORS",
    "SCRFU_LONGITUDINAL_METADATA",
    "SCRFU_LONGITUDINAL_OUTDIR",
    "SCRFU_LONGITUDINAL_SAMPLE_KEY",
    "SCRFU_LONGITUDINAL_DONOR_KEY",
    "SCRFU_LONGITUDINAL_TIME_KEY",
]
missing = [name for name in required if not os.environ.get(name)]
if missing:
    raise SystemExit(f"Unset required variables: {missing}")

def read_table(path):
    sep = "\t" if ".tsv" in path.name.lower() else ","
    return pd.read_csv(path, sep=sep, compression="infer")

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

receptor_path = Path(os.environ["SCRFU_LONGITUDINAL_RECEPTORS"]).expanduser().resolve()
metadata_path = Path(os.environ["SCRFU_LONGITUDINAL_METADATA"]).expanduser().resolve()
outdir = Path(os.environ["SCRFU_LONGITUDINAL_OUTDIR"]).expanduser().resolve()
for label, path in (("receptors", receptor_path), ("metadata", metadata_path)):
    if not path.is_file():
        raise SystemExit(f"Configured {label} path is not a file: {path}")

sample_key = os.environ["SCRFU_LONGITUDINAL_SAMPLE_KEY"]
donor_key = os.environ["SCRFU_LONGITUDINAL_DONOR_KEY"]
time_key = os.environ["SCRFU_LONGITUDINAL_TIME_KEY"]
compartment_key = os.environ.get("SCRFU_LONGITUDINAL_COMPARTMENT_KEY") or None
receptors = read_table(receptor_path)
metadata = read_table(metadata_path)
scrfu.pp.validate_receptor_table(receptors, strict=True)
if sample_key not in receptors:
    raise SystemExit(f"Receptors lack declared sample key: {sample_key}")
design = scrfu.tl.validate_longitudinal_design(
    metadata,
    sample_key=sample_key,
    donor_key=donor_key,
    time_key=time_key,
    compartment_key=compartment_key,
)
receptor_samples = set(receptors[sample_key].dropna().astype(str))
metadata_samples = set(design.design_table[sample_key].dropna().astype(str))
if receptor_samples - metadata_samples:
    raise SystemExit("One or more receptor samples are absent from metadata.")

report = {
    "schema_version": "1",
    "receptor_sha256": sha256(receptor_path),
    "metadata_sha256": sha256(metadata_path),
    "receptor_rows": len(receptors),
    "unique_cdr3": int(receptors["cdr3aa"].dropna().astype(str).nunique()),
    "samples": len(design.design_table),
    "donors": len(design.ordered_donors),
    "timepoints": len(design.ordered_timepoints),
    "compartments": (
        int(design.design_table[compartment_key].nunique()) if compartment_key else None
    ),
    "design_warning_count": len(design.warnings),
    "sample_key": sample_key,
    "donor_key": donor_key,
    "time_key": time_key,
    "compartment_key": compartment_key,
    "scrfu_version": scrfu.__version__,
    "contains_identifiers": False,
}
outdir.mkdir(parents=True, exist_ok=True)
(outdir / "input_validation.json").write_text(
    json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
print(json.dumps(report, sort_keys=True))
PY
```

After this command passes, provide the three required variables, explicit key
names, and the generated `input_validation.json` to the next governed execution
session. That session will run the already frozen RFU, comparator,
longitudinal, and donor-level resampling definitions in
[`methods_freeze.md`](methods_freeze.md). Until then, Figure 2 and its Results
claim remain blocked.
