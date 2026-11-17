# Executable synthetic quickstart

The default tutorial is fully synthetic and redistributable. It validates a
canonical TRB table, uses clearly labelled mock RFU assignments, calculates
repertoire metrics, pseudobulk, phenotype coupling and longitudinal summaries,
matches a tiny synthetic VDJdb-like reference, and exercises experimental BCR
preprocessing and feature extraction. It does not contain RFU assets, VDJdb records, public-cohort rows
or private data.

The frozen fixture SHA256 values are:

- `tutorial_receptors.tsv`:
  `8248399cff8893214b0d4fa7e82f0fcecf4405f20aa1c3641b57cd976f246c09`
- `tutorial_vdjdb.tsv`:
  `42f13868f8f1e9ffc1bc1541c53914149c7af6a4e87c6879bf09a7fad10a0d24`
- `tutorial_bcr.tsv`:
  `8020f24a4afbb81457590b81bb892ada71f966a0c13af9b3779385e76afa2492`

After installing the software candidate, run from any directory:

```bash
python -m scrfu.tutorial --outdir /tmp/scrfu-tutorial
```

The default output manifest states that the assignments are synthetic and must
not be interpreted as official RFU results.

Expected outputs are 12 rows in `assigned_receptors.tsv`, six sample rows in
`rfu_pseudobulk_proportions.tsv`, and a `run_manifest.json` that records
`backend: "mock"`, `external_rfu: false` and `synthetic: true` for every fixture.
The output directory also contains phenotype coupling, longitudinal similarity,
donor retrieval, synthetic VDJdb evidence and experimental BCR feature tables.
These outputs demonstrate software behavior; they are not biological findings.

For real assignment with the external R implementation and model files, follow
[first real RFU run](first_real_rfu.md). The compatibility launcher
`examples/tutorial_end_to_end.py` invokes this same packaged walkthrough.
