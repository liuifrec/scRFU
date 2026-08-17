# Exact Antigen Evidence

`scrfu.tl.annotate_vdjdb()` matches receptor or RFU-result rows to a loaded
`VDJdbReference`. It returns a long evidence table: one query input row by one
reference record. The function never modifies or expands the original RFU
result table. Expansion is isolated in the evidence table and validated against
the expected sequence-to-row reconstruction.

Supported exact modes are:

- `cdr3`: normalized amino-acid CDR3 equality, tier `cdr3_exact`;
- `cdr3_v`: CDR3 plus V-gene equality under the requested V mode, tier
  `cdr3_v_exact`;
- `paired_exact`: exact paired CDR3 information, and paired V when present in
  both inputs, tier `paired_exact`.

Chain-compatible records are required. `chain="TRB"` is the default. Match
tiers describe database matching evidence, not certainty of antigen
specificity. Distance or motif similarity is not labelled exact evidence.

The long output preserves `unique_sequence_id`, `input_row_id`, `cell_id`,
`source_row_id`, RFU assignment fields when present, the query and matched
receptor fields, `reference_row_id`, antigen/MHC/evidence fields, evidence tier,
reference release, and SHA256. Multiple epitopes, MHC alleles, publications, or
database records for one sequence remain separate rows.

`summarize_vdjdb_evidence()` provides derived per-sequence and per-input-row
tables. It reports evidence-record and distinct-epitope counts, maximum score,
observed tiers, and an ambiguity flag. A dominant epitope is populated only
when exactly one distinct epitope is present. The long evidence table remains
authoritative.

`summarize_antigen_context()` joins independent metadata through an explicit,
unique cell identifier. It reports sequence prevalence across samples,
RFU-antigen recurrence across samples/donors, phenotype abundance, evidence
tiers, ambiguity, and join QC. These are descriptive summaries; cells are not
treated as independent biological replicates and no clinical association is
inferred.

## Offline workflow

`examples/vdjdb_antigen_evidence.py` accepts RFU per-sequence results, optional
per-row results and metadata, a local reference, explicit release/checksum, and
all matching/coherence/null policies. It writes reference QC, authoritative
long matches, sequence/row summaries, RFU and global coherence, grouping
baselines, a permutation summary (and optionally values), context tables,
generic figures, and a checksummed run manifest. The workflow has no network
path and never modifies the reference.

```bash
python examples/vdjdb_antigen_evidence.py \
  --rfu-sequences results/rfu_results_per_sequence.tsv.gz \
  --vdjdb /path/to/pinned-vdjdb.tsv.gz \
  --vdjdb-release EXPLICIT_RELEASE_LABEL \
  --expected-sha256 EXPECTED_SHA256 \
  --outdir results/vdjdb_evidence
```
