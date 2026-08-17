# Journal Positioning

## Why This Is More Than a Wrapper

scRFU bridges established RFU-based TCR functional annotation into
single-cell AnnData workflows. The contribution is not a new RFU algorithm; it
is an integration framework that makes sequence-derived functional units usable
alongside transcriptomic metadata, sample annotations, cell states, and
downstream scverse analysis.

The package provides dataset-independent receptor adapters, expression-free
caches, exact-sequence deduplication, restartable execution, explicit
nearest/threshold assignment, convergence and phenotype coupling,
biological-sample pseudobulk, overlap, conventional repertoire comparisons,
visualization, and provenance. These components turn RFU
assignment from a standalone repertoire step into an AnnData-compatible analysis
layer that can be benchmarked on public single-cell immune profiling datasets.

The version-pinned antigen-evidence layer adds a falsifiable validation
question: whether exact VDJdb-labelled sequences are more label-coherent within
RFUs than in size-preserving or receptor-property-stratified null groupings. It
preserves ambiguous database evidence and compares RFU grouping with simple
TRBV and CDR3-length baselines. This is external evidence validation, not a new
antigen-specificity predictor.

## Additional Evidence Needed for Higher-End Submission

- At least two real datasets with completed RFU runs and saved manifests.
- Comparison to conventional repertoire metrics such as clonotype abundance,
  diversity, V/J usage, and cell-type composition.
- Runtime and resource summaries on synthetic and real inputs.
- Robust user documentation and failure-mode guidance.
- Stable v0.x release with fixed package metadata and tagged manuscript code.
- Reproducible figure scripts that rebuild manuscript outputs from result
  directories.
- Clear licensing and attribution statement for the upstream RFU dependency.

## Candidate Journal Tiers

These are possible fits rather than predictions.

- Bioinformatics: plausible if real-data demonstrations, API stability, and
  comparison analyses are strong.
- Briefings in Bioinformatics: possible only with a broader methods/software
  narrative and substantial biological benchmarking.
- NAR Genomics and Bioinformatics: plausible for a reproducible software
  application note or methods-focused submission.
- GigaScience: possible if workflows, data provenance, reproducibility, and
  reusable outputs are emphasized.
- Genome Biology Methods-style angle: only realistic if the biological
  demonstration becomes strong and broadly useful.
- Current Bioinformatics or Computational Biology and Chemistry: fallback
  options if the contribution remains primarily a lightweight integration
  package.

## Risks

- Reviewers may view the package as only a wrapper around upstream RFU.
- Upstream RFU dependency and licensing may complicate packaging or
  reproducibility.
- Lack of a novel algorithm may reduce appeal for higher-tier venues.
- Biological impact will be weak without completed real-data demonstrations.
- Incomplete support for all scirpy AIRR storage variants may be flagged.

## Mitigation

- Emphasize AnnData/scirpy integration, reproducible workflows, provenance, and
  downstream scverse compatibility.
- Include at least two real datasets and transparent benchmark manifests.
- Compare RFU outputs with clonotype/diversity/V-gene metrics rather than
  presenting RFU matrices in isolation.
- Keep upstream RFU attribution and dependency boundaries explicit.
- Maintain CI-safe tests, small examples, stable API documentation, and
  manuscript-ready output schemas.

## Explicit Future Work

Broader antigen-coherence validation across independently selected datasets,
exploratory sequence-similarity candidates, RFU networks, modeled differential
abundance, BCR lineage utilities, and experimental ALFU remain future
directions. Current outputs make no definitive antigen-specific, clinical,
radiation, or B-cell functional claims.
