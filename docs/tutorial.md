# End-to-end public tutorial

The default tutorial is fully synthetic and redistributable. It validates a
canonical TRB table, uses clearly labelled mock RFU assignments, calculates
repertoire metrics, pseudobulk, phenotype coupling and longitudinal summaries,
matches a tiny synthetic VDJdb-like reference, and exercises experimental BCR
preprocessing. It does not contain RFU assets, VDJdb records, public-cohort rows
or private data.

The frozen fixture SHA256 values are:

- `tutorial_receptors.tsv`:
  `8248399cff8893214b0d4fa7e82f0fcecf4405f20aa1c3641b57cd976f246c09`
- `tutorial_vdjdb.tsv`:
  `42f13868f8f1e9ffc1bc1541c53914149c7af6a4e87c6879bf09a7fad10a0d24`
- `tutorial_bcr.tsv`:
  `8020f24a4afbb81457590b81bb892ada71f966a0c13af9b3779385e76afa2492`

From a source checkout with scRFU installed:

```bash
python examples/tutorial_end_to_end.py --outdir /tmp/scrfu-tutorial
```

The default output manifest states that the assignments are synthetic and must
not be interpreted as official RFU results.

To exercise the real canonical TRB backend, supply an untouched external RFU
checkout explicitly:

```bash
python examples/tutorial_end_to_end.py \
  --backend rfu_repo \
  --rfu-dir /path/to/external/RFU \
  --outdir /tmp/scrfu-tutorial-official
```

The real-backend command requires R and the upstream assets; they are never
bundled with scRFU. BCR processing in this tutorial is experimental and never
invokes the TCR RFU backend.
