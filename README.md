# scRFU

scRFU is a transferable, frozen-reference receptor-state framework for
longitudinal, cross-cohort, and single-cell repertoire analysis. It preserves
the established RFU assignment method while adding dataset-independent receptor
interfaces, deterministic execution, cell-level reconstruction, reference-
coverage diagnostics, and sample-level analysis contracts. Wells is a named
adapter and bounded public benchmark, not the data model.

## What scRFU Does

- Reads current Scirpy Awkward AIRR records directly from AnnData or a MuData
  AIRR modality while preserving every receptor chain and its order.
- Prepares canonical receptor rows through Wells, scirpy/AIRR, generic
  DataFrame, and Cell Ranger VDJ adapters while keeping metadata separate.
- Selectively reads receptor and observation data without loading expression
  matrices and writes portable checksummed receptor caches.
- Calls a user-provided upstream RFU repository backend.
- Attaches RFU labels, RFU scores, extracted TRB features, and provenance into
  AnnData.
- Provides RFU summaries, group-level RFU matrices, matplotlib plots, and TSV
  export for downstream analysis.
- Provides conventional repertoire metrics, sample-level RFU pseudobulk,
  overlap, phenotype coupling, explicit convergence definitions, and generic
  manuscript-oriented plots.
- Provides CI-safe synthetic examples plus real-data workflow scaffolds for
  user-provided datasets.
- Triangulates user-supplied RFU-QTL, eQTL, caQTL, and GWAS tables with build-aware
  variant matching, descriptive credible-set overlap, and explicit provenance.
  See [regulatory triangulation](docs/regulatory_triangulation.md) for the
  experimental API and local-file Matos CD4 example. Overlap does not establish
  statistical colocalization, causality, or antigen specificity.
- Matches RFU sequences to user-supplied, version-labelled local VDJdb tables
  and benchmarks descriptive antigen-label coherence against explicit nulls.
- Validates repeated-measures designs and constructs longitudinal RFU matrices,
  pairwise similarity tables, donor-retrieval results, and explicit descriptive
  trajectory classes without imputing unobserved visits.
- Records immutable frozen-reference manifests and produces target-cohort
  coverage, pseudobulk, and harmonization diagnostics without refitting the RFU
  reference in the target cohort.
- Provides explicitly experimental BCR canonicalization, heavy/light pairing,
  QC, and conservative receptor-state features without applying TCR RFUs.

## What scRFU Does Not Do

- scRFU does not invent RFU or replace the upstream RFU method.
- scRFU does not vendor upstream RFU code, atlas files, or `.Rdata` objects.
- scRFU does not download real datasets.
- scRFU does not bundle or silently download VDJdb, and a database match is not
  proof that an RFU is antigen-specific.
- Experimental BCR preprocessing is not a validated BCR functional-unit model;
  no frozen BCR reference is shipped or claimed.
- Current longitudinal and transfer APIs are methodological foundations; they
  are not population-level evidence until evaluated on independent cohorts.
- Accepting IGH/IGK/IGL rows in the receptor schema does not constitute a BCR
  functional-unit model.

## Installation

Install from a source checkout:

```bash
git clone https://github.com/liuifrec/scRFU.git
cd scRFU
python -m pip install .
```

Or install a reviewed local wheel:

```bash
python -m pip install /path/to/scrfu-0.4.0-py3-none-any.whl
```

The project does not document `pip install scrfu` as an available distribution
until a PyPI upload has actually occurred. For development:

```bash
python -m pip install -e ".[dev]"
```

Optional scirpy extras:

```bash
python -m pip install -e ".[dev,scirpy]"
```

Core imports do not require matplotlib. Install plotting support with
`python -m pip install ".[plotting]"`; `.[all]` installs plotting, MuData, and
scirpy extras.

## External RFU Dependency

RFU assignment requires a separate checkout of the upstream RFU repository.
Pass that location with `rfu_dir`/`--rfu-dir`, or set `RFU_DIR`. An explicit
argument takes precedence over the environment variable. Standard mode uses the
official public `AssignRFUs()` implementation and remains the default. scRFU
records hashes and detected optional capabilities, but the files remain outside
this repository.

Use `scrfu doctor` to inspect Python, optional dependencies, R availability,
external RFU configuration, artifact hashes, and output writability. Paths are
shortened by default; `scrfu doctor --verbose` prints full runtime paths.

## Quickstart

For a runnable data-free walkthrough, see
[`docs/tutorial.md`](docs/tutorial.md) and the
[scverse-native tutorial](docs/scverse_tutorial.md). The default tutorials use synthetic
fixtures and clearly labelled mock assignments; real RFU execution remains an
explicit external-backend step.

```python
import scrfu

# Primary current Scirpy/AnnData/MuData path. Results are aligned to every
# source AIRR chain in obsm["scrfu"]. Cell-level obs values are opt-in.
scrfu.tl.assign_rfu(
    mdata,
    rfu_dir="/path/to/RFU",
    airr_mod="airr",
    airr_key="airr",
    key_added="scrfu",
    cell_summary=True,
    summary_policy="ambiguity_aware",
)

prepared = scrfu.adapters.prepare_receptors(
    airr_dataframe,
    adapter="generic_airr_dataframe",
    chain="TRB",
)
table_run = scrfu.tl.call_rfu_table(
    prepared.receptors,
    rfu_dir="/path/to/RFU",
    chunk_size=20_000,
    max_workers=2,
)

qc = scrfu.tl.validate_airr(adata, airr_key="airr", chain="TRB")

# Historical one-row-per-cell compatibility path:
scrfu.tl.call_rfu(
    adata,
    backend="rfu_repo",
    rfu_dir="/path/to/RFU",
    mode="standard",
    chunk_size=20_000,
    workdir="results/rfu_work",
)

summary = scrfu.tl.rfu_summary(adata, groupby="sample")
matrix = scrfu.tl.aggregate_rfu(adata, groupby="sample")

pseudobulk = scrfu.tl.rfu_pseudobulk(
    per_row_results,
    sample_key="donor_id",
    weighting="cell",
    normalize="proportion",
)
overlap = scrfu.tl.rfu_overlap(pseudobulk, metric="jaccard")

design = scrfu.tl.validate_longitudinal_design(
    per_row_results,
    sample_key="sample_id",
    donor_key="donor_id",
    time_key="timepoint",
)
longitudinal = scrfu.tl.rfu_longitudinal_matrix(
    per_row_results,
    sample_key="sample_id",
    donor_key="donor_id",
    time_key="timepoint",
    design=design,
)
pairs = scrfu.tl.longitudinal_similarity(longitudinal, metric="cosine")

reference = scrfu.tl.load_vdjdb_reference(
    "vdjdb-release.tsv.gz",
    release_label="EXPLICIT_RELEASE_LABEL",
    expected_sha256="EXPECTED_SHA256",
)
evidence = scrfu.tl.annotate_vdjdb(per_sequence_results, reference)
coherence = scrfu.tl.rfu_antigen_coherence(per_sequence_results, evidence)

scrfu.pl.rfu_bar(adata, groupby="sample")
scrfu.pl.rfu_heatmap(adata, groupby="sample")
scrfu.pl.rfu_score_hist(adata)
```

## Native AIRR storage

The primary path reads current Scirpy records from `adata.obsm["airr"]` or
`mdata.mod["airr"].obsm["airr"]`. Chain-level results are written to an aligned
Awkward array in `obsm["scrfu"]`; RFU provenance and schema version are stored
in `uns["scrfu"]` of the AIRR modality. TRA, BCR, nonproductive, malformed, and
missing-sequence chains remain explicit non-eligible records. Multiple TRBs are
never collapsed. H5AD/H5MU round trips and observation subsetting preserve the
alignment. Use `scrfu.tl.validate_scrfu_schema` after loading and
`scrfu.tl.concat_scrfu` when concatenating annotated AnnData objects.

Cell-level RFU columns are created only when requested. The default summary is
ambiguity-aware and leaves the cell RFU missing when eligible TRB chains map to
different RFUs. `scirpy.pp.index_chains` is not required for chain assignment;
it is used only by the explicit `primary_vdj` summary policy.

The compatibility table path accepts AIRR-like aliases including
cell barcode (`cell_id`, `cell`, `barcode`, `cellid`), chain (`chain`, `locus`),
CDR3 amino-acid sequence (`cdr3aa`, `junction_aa`, `cdr3_aa`, `cdr3`), V gene
(`v_call`, `v_gene`, `trbv`, `v`), and optional `productive`.

The historical `call_rfu` wrapper writes one-row-per-cell compatibility fields:

- `adata.obs["trb_cdr3aa"]`
- `adata.obs["trbv"]`
- `adata.obs["rfu_label"]`
- `adata.obs["rfu_score"]`
- `adata.uns["scrfu"]`

See [docs/scverse_storage_schema.md](docs/scverse_storage_schema.md) for the
native schema and [docs/api_contract.md](docs/api_contract.md) for compatibility
contracts.

The canonical schema, adapters, and expression-free cache are documented in
[docs/receptor_schema.md](docs/receptor_schema.md),
[docs/adapters.md](docs/adapters.md), and
[docs/receptor_cache.md](docs/receptor_cache.md). Expression matrices are not
required after receptor preparation.

Reference coverage is the fraction of eligible receptor rows whose nearest RFU
score passes the configured frozen-reference threshold. Nearest assignment and
threshold-qualified assignment are separate policies: a threshold miss is a
low-reference-similarity result, not a calibrated confidence probability or an
automatic out-of-distribution declaration.

## Analysis Utilities

The public analysis layer includes legacy summaries plus
`repertoire_metrics`, `rfu_metrics`, `rfu_pseudobulk`, `rfu_overlap`,
`rfu_phenotype_coupling`, and generic `scrfu.pl` functions. See
[docs/v0.2_api_review.md](docs/v0.2_api_review.md) for the compatibility
boundary and the focused metric documentation for formulas.

The offline antigen-evidence API is documented in
[docs/vdjdb_reference.md](docs/vdjdb_reference.md),
[docs/antigen_evidence.md](docs/antigen_evidence.md),
[docs/rfu_antigen_coherence.md](docs/rfu_antigen_coherence.md), and
[docs/antigen_null_models.md](docs/antigen_null_models.md). Exact matching tiers
describe evidence strength, not biological certainty.

## Reproducible Examples

- [examples/synthetic_scirpy_demo.py](examples/synthetic_scirpy_demo.py):
  offline synthetic AnnData/scirpy-style demo.
- [examples/tutorial_scverse_native.py](examples/tutorial_scverse_native.py):
  current Scirpy/MuData chain-aligned tutorial with explicitly synthetic mock
  assignments.
- [examples/real_scirpy_workflow.py](examples/real_scirpy_workflow.py):
  template for user-provided real h5ad input.
- [examples/wells_atlas_workflow.py](examples/wells_atlas_workflow.py):
  restartable Wells public-atlas workflow for a user-supplied H5AD or compact
  expression-free receptor cache.
- [examples/receptor_table_workflow.py](examples/receptor_table_workflow.py):
  dataset-independent cache, canonical TSV, or named-H5AD-adapter workflow.
- [examples/rfu_downstream_analysis.py](examples/rfu_downstream_analysis.py):
  sample-level downstream metrics and generic figures from RFU per-row results.
- [examples/benchmark_summary.py](examples/benchmark_summary.py): combines
  user-labelled manifests into a manuscript-ready benchmark table.
- [examples/vdjdb_antigen_evidence.py](examples/vdjdb_antigen_evidence.py):
  offline exact antigen evidence, ambiguity-aware coherence, simple baselines,
  permutation null, and generic figures from a pinned local reference.

For large inputs, chunking occurs after exact-CDR3 deduplication. Completed
chunks are reused only after their manifest, scientific inputs, hashes, output
schema, row count, and identifiers all validate. `chunk_size=None` retains the
existing one-call behavior; `resume=False` reruns every chunk, and
`force_recompute=True` takes precedence and forces every chunk to run again.
The default remains serial. With explicit `max_workers > 1`, independent chunks
run concurrently in isolated directories and final results are concatenated in
chunk-index order. Worker count and executor are recorded in provenance.

## Manuscript and Benchmark Workflows

- [examples/benchmark_scirpy_dataset.py](examples/benchmark_scirpy_dataset.py):
  benchmark scaffold for user-provided AnnData/scirpy h5ad files.
- [examples/gse190905_prepare.py](examples/gse190905_prepare.py):
  preparation scaffold for user-downloaded GSE190905-like TCR and metadata
  files.
- [docs/benchmark_outputs.md](docs/benchmark_outputs.md): output schema.
- [docs/repository_audit.md](docs/repository_audit.md): manuscript-readiness
  audit.
- [docs/manuscript_figures.md](docs/manuscript_figures.md): figure-generation
  plan.

## Citation and Attribution

Please cite and attribute the original upstream RFU method when using RFU
assignments. scRFU is an integration framework around user-provided upstream RFU
code/data. Cite scRFU separately using the metadata in
[`CITATION.cff`](CITATION.cff); no DOI or publication citation is claimed here.

## Documentation and License

The runnable tutorial is in [`docs/tutorial.md`](docs/tutorial.md), the API
documentation sources are under [`docs/api/`](docs/api/), and configuration is
documented in [`docs/configuration.md`](docs/configuration.md). scRFU's own code
is distributed under the [MIT license](LICENSE). External RFU assets, VDJdb,
Wells, GEO, OAS, and other public datasets retain their separate source terms
and are not redistributed by the package.
