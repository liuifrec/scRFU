# Month 2 execution report

Execution date: 2026-08-24. Tracked code contains no runtime absolute paths or
result tables. Runtime artifacts are external under `SCRFU_MONTH2_OUTDIR`; this
report uses paths relative to that variable. No commit, tag, release, or push
was made.

1. **Initial test status.** Branch `feature/month2-real-validation`; initial
   worktree clean; Python 3.10.20; scRFU 0.1.0; 261 tests passed and 4 skipped;
   Ruff lint/format and `git diff --check` passed.
2. **Package metadata fixes.** Author, email, homepage, repository, and issues
   metadata were already correct, with no placeholder strings. Version remains
   0.1.0. The documented coherent next release remains 0.4.0 because historical
   v0.2.0/v0.3.0 tags exist.
3. **Official RFU checkout/version/hash.** Untouched public upstream commit
   `ff7a1ea6aca28444a7ff8ec7ef09745e7341f99d`; `RFU.R` SHA256
   `92c2faa33f2e7f60d6470ad7dfd653eb0ca54859aa198ed53ad0c876aae8640b`
   (20,161 bytes), trimer reference
   `820fb71428913974e994543cbbbaa591f54355cdb73d651f068f0e691a6cdffd`
   (7,292,940 bytes), centroid reference
   `64783074360edeca84b3f49ab9d682263e954752fc09bb208edaabde0fd82553`
   (28,253,173 bytes). Source: `https://github.com/s175573/RFU`; recorded
   2026-08-24. No asset was copied into scRFU.
4. **Real original-RFU parity.** Twelve shuffled adversarial rows, five unique
   eligible CDR3s, ten eligible rows, and two noneligible rows. RFU ID 10/10,
   labels, scores, threshold status, reconstruction, and order all matched;
   maximum absolute score difference 0 at tolerance `1e-12`. Total backend time
   3.924 s; child peak RSS 352,116 KB. Evidence:
   `original_rfu_parity/parity_summary.json`.
5. **Wells subset construction.** Selective HDF reads only; no `X` or `raw/X`.
   Fixed seed 20260824 and deterministic random sampling preserved source order.
   Source: 610,429 cells, 24 donors, 17 cell types, 10 tissues, 247 libraries.
   A proportional donor×tissue×cell-type 25k check reduced stratum total
   variation from 0.0772 to 0.0115; neither sample is claimed representative.
6. **1k/10k/25k receptor counts.** Respectively 499/4,861/12,415 productive
   primary-TRB rows and 453/4,130/10,038 unique CDR3 queries. Selected-cell hashes
   and full metadata counts are in `wells_bounded/wells_sampling_summary.tsv`.
7. **Serial scaling.** Unchunked serial wall times were 3.920 s, 18.141 s, and
   40.276 s for 1k, 10k, and 25k subsets.
8. **Parallel scaling.** At 25k, 5,000-query chunks took 45.650 s serial and
   25.859 s with two process workers; 1,000-query chunks took 63.809 s and
   36.102 s. Four workers were not required.
9. **Restart/cache benefit.** At 25k, resumed runs were 0.210–0.222 s versus
   25.859–63.809 s fresh. Held-out GSE157007 reused all nine chunks in 0.871 s
   after a 117.300 s fresh run.
10. **Peak-memory results.** Wells 25k Python peak RSS was 482,296 KB and the
    largest reported RFU backend peak was 1,802,608 KB. GSE157007 fresh backend
    peak was 947,596 KB.
11. **Deduplication benefit.** Unique-query/input ratios were 0.908, 0.850, and
    0.809 for Wells 1k/10k/25k. GSE190905 and GSE157007 ratios were 0.593 and
    0.717.
12. **Chunk-size invariance.** Every meaningful Wells configuration produced
    the same assignment SHA256; 20k chunks were skipped because each subset had
    fewer than 20k unique queries.
13. **Input-order invariance.** The shuffled Wells 25k run exactly matched the
    source-order result.
14. **Downsampling robustness.** For nearest-policy cell subsampling, median
    cosine similarity to the 25k reference was 0.526/0.730/0.885/1.000 at
    25/50/75/100%; corresponding median Spearman values were approximately
    0.510/0.712/0.867/1.000. Full metrics, seeds, undefined constant-vector
    correlations, sequence subsampling, and multinomial resampling are retained
    in `wells_bounded/figure1_robustness_source.tsv`.
15. **Reference coverage.** Wells 25k threshold-pass fraction was 0.7694;
    GSE190905 0.7831; held-out GSE157007 0.7688. Nearest and threshold-qualified
    policies are reported separately.
16. **Wells phenotype coupling.** Nearest-policy cell-subsample coupling cosine
    means were 0.506/0.708/0.863/1.000 at 25/50/75/100%, with dominant-phenotype
    agreement 0.637/0.735/0.853/1.000. Threshold-policy results were similar.
17. **Conventional repertoire comparisons.** Wells completed exact CDR3, V,
    CDR3 length, Shannon, and Simpson; J/clonotype were absent. GSE190905 and
    GSE157007 completed exact CDR3, clonotype, V, J, length, Shannon, and
    Simpson. Edit-distance clustering was skipped when unique CDR3 count
    exceeded the frozen quadratic limit of 2,000.
18. **VDJdb real-reference status/results.** Blocked: `VDJDB_PATH` and
    `VDJDB_RELEASE` are unset. Nothing was downloaded. The exact deferred run
    commands are below.
19. **Longitudinal private-data input status.** All
    `SCRFU_LONGITUDINAL_*` variables are unset. No filesystem search occurred.
20. **Longitudinal results if configured.** Not run. The six-volunteer cohort
    remains a deep methods-validation cohort, not population evidence.
21. **Public dataset candidates.** PRJNA602091 is the development aging cohort;
    GSE190905 the first independent technical validation; GSE157007 the held-out
    aging/frailty cohort; GSE158848 the reserve. See
    [`public_dataset_candidates.md`](public_dataset_candidates.md).
22. **Public datasets acquired.** GSE190905 receptor/metadata files and 17
    receptor-only GSE157007 files (4.5 MB total) were acquired externally and
    recorded by hash. Large expression archives and PRJNA602091 raw reads were
    not downloaded.
23. **Cross-cohort transfer results.** GSE190905 contained 27,655 receptor rows,
    16,391 unique CDR3s, 4,465 nearest RFUs, and 3,940 threshold-qualified RFUs.
    Leave-one-timepoint-out RFU cosine retrieval was top-1 0.667, top-3 0.833,
    MRR 0.764. Held-out GSE157007 contained 60,125 receptor rows, 43,082 unique
    CDR3s, 4,898 nearest RFUs, and 4,618 threshold-qualified RFUs.
24. **Held-out validation status.** Completed for applicable prespecified
    technical endpoints using the immutable input hash
    `6e311d8c537da80f1aa023e06c5708c58ebbf4e4ce902415d98074db1ca75732`.
    RFU nearest cosine stability averaged 0.936 at 50% and 0.976 at 75% cell
    subsampling. Donor retrieval was undefined because each donor has one sample.
25. **Comparator results.** GSE190905 used identical leave-one-timepoint-out
    partitions. RFU, exact CDR3, and clonotype cosine top-1 were each 0.667;
    their top-3 values were 0.833, 0.750, and 0.833. The result supports
    complementary trade-offs, not universal RFU superiority. GSE157007
    comparator stability is in its held-out source tables; sample-local
    clonotype IDs were namespaced before cross-sample comparison.
26. **Figure 1 source-table status.** Complete for bounded parity, scaling,
    memory, resume, invariance, deduplication, threshold, and robustness.
27. **Figure 2 source-table status.** Blocked because explicit longitudinal
    inputs are not configured.
28. **Figure 3 source-table status.** Partial: Wells phenotype coupling,
    GSE190905 transfer/comparators, and GSE157007 held-out transfer/stability are
    present; VDJdb and PRJNA602091 aging development evidence are absent.
29. **Manuscript evidence-matrix status.** Updated in
    [`manuscript_evidence_matrix.md`](manuscript_evidence_matrix.md), separating
    executed evidence from infrastructure and blockers.
30. **Submission-gate changes.** Original parity, runtime/memory, and held-out
    technical validation are complete; cross-cohort/comparator/source-table
    gates remain partial; longitudinal, VDJdb, release, and DOI remain blocked.
    BCR is removed from the first manuscript.
31. **Tests/Ruff/build results.** Final suite: 262 passed, 4 skipped. Official
    RFU integration: 3 passed, 1 capability-dependent skip; secondary modified
    checkout: 3 passed, 1 skip. `ruff check .`, `ruff format --check .`, and
    `git diff --check` passed. Isolated sdist/wheel build passed after network
    access was allowed for the declared hatchling dependency. Wheel inspection
    found no RFU assets, H5AD, VDJdb, Rdata, private paths, or credentials.
32. **Files modified.** Dataset/held-out/submission-gate documentation; three
    acquisition manifests; public-candidate, evidence-matrix, and this report;
    original-RFU parity, Wells Month 2, GSE190905 preparation, and GSE157007
    held-out examples; and the GSE190905 preparation tests. All data and result
    tables remain outside the repository.
33. **Scientific blockers.** No real VDJdb reference; no configured governed
    longitudinal cohort; PRJNA602091 preprocessing/acquisition unresolved; no
    population-level interpretation is justified from GSE190905 or bounded Wells.
34. **Software blockers.** No remaining demonstrated execution blocker. The
    held-out run exposed and fixed example-only adapter and wrapper-path errors;
    Cell Ranger sample-local clonotype identifiers required explicit namespacing.
35. **Exact next execution step.** Configure the version-pinned VDJdb variables
    and run the commands below; if unavailable, freeze the PRJNA602091 run/size
    manifest and UMI-aware preprocessing workflow before downloading raw reads.

## Deferred VDJdb command

```bash
export SCRFU_MONTH2_OUTDIR=/path/to/month2-output
python -c 'import os, pandas as pd; root=os.environ["SCRFU_MONTH2_OUTDIR"]; p=f"{root}/wells_bounded/subsets/wells_25000/rfu_results_per_row.tsv.gz"; d=pd.read_csv(p,sep="\t"); d.drop_duplicates("unique_sequence_id").to_csv(f"{root}/wells_bounded/wells_25000_unique_sequences.tsv.gz",sep="\t",index=False)'
for match in cdr3 cdr3_v; do
  for policy in nearest threshold_pass; do
    for ambiguity in fractional exclude_ambiguous; do
      python examples/vdjdb_antigen_evidence.py \
        --rfu-sequences "$SCRFU_MONTH2_OUTDIR/wells_bounded/wells_25000_unique_sequences.tsv.gz" \
        --rfu-rows "$SCRFU_MONTH2_OUTDIR/wells_bounded/subsets/wells_25000/rfu_results_per_row.tsv.gz" \
        --vdjdb "$VDJDB_PATH" --vdjdb-release "$VDJDB_RELEASE" \
        --match-mode "$match" --assignment-policy "$policy" \
        --ambiguity-policy "$ambiguity" --n-permutations 1000 \
        --random-state 20260824 \
        --outdir "$SCRFU_MONTH2_OUTDIR/vdjdb/${match}_${policy}_${ambiguity}"
    done
  done
done
```

## Deferred private-input validation command

Set all three paths and explicitly set the key names; no default is silently
used when the schema is ambiguous.

```bash
export SCRFU_LONGITUDINAL_RECEPTORS=/path/to/receptors.tsv.gz
export SCRFU_LONGITUDINAL_METADATA=/path/to/metadata.tsv.gz
export SCRFU_LONGITUDINAL_OUTDIR=/path/to/private-output
export SCRFU_LONGITUDINAL_SAMPLE_KEY=sample_id
export SCRFU_LONGITUDINAL_DONOR_KEY=donor_id
export SCRFU_LONGITUDINAL_TIME_KEY=time
export SCRFU_LONGITUDINAL_COMPARTMENT_KEY=compartment
python - <<'PY'
import os
import pandas as pd
import scrfu

def read(path):
    return pd.read_csv(path, sep="\t" if ".tsv" in path else ",")

required = [
    "SCRFU_LONGITUDINAL_RECEPTORS", "SCRFU_LONGITUDINAL_METADATA",
    "SCRFU_LONGITUDINAL_OUTDIR", "SCRFU_LONGITUDINAL_SAMPLE_KEY",
    "SCRFU_LONGITUDINAL_DONOR_KEY", "SCRFU_LONGITUDINAL_TIME_KEY",
    "SCRFU_LONGITUDINAL_COMPARTMENT_KEY",
]
missing = [name for name in required if not os.environ.get(name)]
if missing:
    raise SystemExit(f"Unset required variables: {missing}")
receptors = read(os.environ["SCRFU_LONGITUDINAL_RECEPTORS"])
metadata = read(os.environ["SCRFU_LONGITUDINAL_METADATA"])
scrfu.pp.validate_receptor_table(receptors, strict=True)
design = scrfu.tl.validate_longitudinal_design(
    metadata,
    sample_key=os.environ["SCRFU_LONGITUDINAL_SAMPLE_KEY"],
    donor_key=os.environ["SCRFU_LONGITUDINAL_DONOR_KEY"],
    time_key=os.environ["SCRFU_LONGITUDINAL_TIME_KEY"],
    compartment_key=os.environ["SCRFU_LONGITUDINAL_COMPARTMENT_KEY"],
)
print({"receptor_rows": len(receptors), "samples": len(design.design_table),
       "donors": len(design.ordered_donors),
       "timepoints": len(design.ordered_timepoints)})
PY
```
