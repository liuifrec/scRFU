# Restored RP1-14 reuse audit

The authorized published-only RP1-14 material was restored locally on
2026-09-16. Both biological input directories are untracked and protected by
exact-path ignore rules; `data/schema.md` remains tracked. Neither sequences,
original sample names nor detailed biological results are redistributed.

## Assets and design

| Restored asset | Files / schema | Interpretation |
| --- | --- | --- |
| Source repertoire exports | Six ZIPs containing 36 tab-separated tables, 46 columns each | 4,342,203 nucleotide-sequence rows; 70,960,689 sequencing reads. Half use `count`, half `count (reads)`. |
| `TCR variable definition.xlsx` | Variable/definition workbook | Documents full nucleotide sequence, CDR3 amino acids, CDR3 nucleotide length, zero-based conserved-cysteine start, V/J calls and read percentages. |
| `TCR-seq gender and age.xlsx` | 36 sample metadata rows | Six donors, three visits each, both CD4 and CD8 at every visit; no missing samples in this restored subset. |
| Historical receptor RFU exports | 36 `*.cdr3_rfu.tsv` tables | `cdr3_aa`, `trbv`, `freq`, `rfu`, `max_cor`, `pass_thr`; 359,530 rows. All saved `freq` values are missing. |
| `RFU_matrix.tsv` | 5,000 rows × 36 sample columns | Unweighted saved assignment rows per 10,000, including below-threshold nearest labels. No logarithm or read weighting. |

Raw files total 155,430,399 compressed/workbook bytes; historical RFU files
total 23,662,403 bytes. Detailed filenames, ZIP members, source sample crosswalks
and file hashes remain external. Source percentages agree with the complete
sample read denominator. Counts are sequencing reads, not independent cells
or experimentally observed templates. Estimated-genome fields are not used as
measured template counts.

Visits are ordered by the source collection-age workbook. File suffixes 3, 2, 1
correspond to chronological visits 1, 2, 3; sorting suffixes would reverse time.
First-to-last follow-up ranges from **14.95 to 24.86 years**, rather than an
assumed identical 20-year interval. Collection ages span 23.50–65.04 years.
CD4 and CD8 collection ages agree within every donor/visit. Aliases RP01–RP06
are used only in external derived tables and figures.

The corresponding prior publication is
[Yoshida et al., Experimental Gerontology (2017)](https://doi.org/10.1016/j.exger.2017.05.015).
It describes six healthy volunteers, not an atomic-bomb survivor exposure
comparison. Its reported long-term persistence and CD4/CD8 differences are
prior knowledge, not newly prespecified aging biomarkers for this analysis.

## Demonstrated export-alignment error and minimal repair

In every sample, the saved sequence/V labels equal the **first N** rows of the
top 10,000 nonmissing translated source sequences. However, the saved assignment
vector has the length of the subset **starting with C**. The upstream
`EncodeRepertoire` function applies that C-start filter. Consequently, saved
labels after a removed sequence were shifted relative to their assignments.
The unrepaired export contains 41,572 repeated amino-acid sequences with
conflicting RFUs, including within-sample conflicts.

Reattaching the existing vector to the corresponding filtered source rows
repairs 330,433 amino-acid row labels. **No RFU number, score or qualification
value is changed.** Exact source-prefix and count checks pass in all 36 files;
every repeated amino-acid sequence then has one RFU and one score. The saved
matrix is reconstructed from those unchanged vectors with maximum absolute
error 4.98 × 10⁻¹⁴. Original biological files remain untouched.

Historical flags agree with correlation threshold 0.6. RFU numbers lie in
1–5,000, with five unassigned rows retained in the audit. Translation of the
source nucleotide CDR3 interval agrees with its supplied amino acids. Main
analyses retain productive, canonical amino-acid sequences; historical stop-
containing translated rows remain in the reconciliation but not the endpoints.

## Reference compatibility and fixed observation universe

The locally recovered reference file matches the current `km5000noMax.Rdata`
SHA256 `64783074360edeca84b3f49ab9d682263e954752fc09bb208edaabde0fd82553`.
A seeded, bounded comparison of six receptors per sample gives **216/216**
identical RFU labels with the current standard backend; maximum score difference
is 1.12 × 10⁻¹⁵. The recovered `RFU.R` contains `AssignRFUs_with_map` and `RFUbatch_with_maps`,
including the exact prefix truncation that explains the export error. Its SHA256
is `47e3d529d1f9c9b7a0370467be1850ba13a1f1117d8da14494b1d2ebbf04a000`.
The original invocation/notebook and contemporaneous execution manifest are
unavailable. This is strong bounded computational compatibility,
not proof that every historical run used an independently verified file hash.
No full assignment campaign was performed.

The repaired amino-acid-to-RFU dictionary is fixed across all samples and applied
to full productive source counts. This reuses assignments when a sequence falls
below another visit's historical top-10,000 cutoff. It does not assign previously
unmapped sequences. Excluded mass is reported separately for every sample;
the main estimand remains conditional on this recoverable reference subset.
Observed receptor identity is translated/verified nucleotide CDR3 plus the
source maximum-resolved V and J calls. It is not a verified paired-chain lineage.

RFU, V and V/J distances use the same receptors and normalization. Family-only
or unresolved gene calls remain explicitly marked categories. Read weighting,
one-per-observed-receptor weighting, removal of dominant receptors and empirical
read subsampling are separate analyses. The original development seeds, 30 fixed
grouping controls, 50 subsamples and depth cap of 500 are retained, with the
sampling unit explicitly changed from single cells to source reads. No calibrated
null p-value or donor confidence interval is inferred from these replicates.

## External records and commands

Preparation lives in
`$SCRFU_METHODS_DIR/prepared/rp1_14_v1/`, including
`rp1_14_asset_manifest.tsv`, `rp1_14_provenance.json`,
`historical_alignment_repair.tsv`, `sample_registry.tsv`,
`source_sample_crosswalk.tsv`, the receptor parquet and 36 source-count parquets.
The bounded compatibility comparison is saved separately. Preparation has a
validated stage checkpoint and a final checksum manifest.

```bash
SCRFU_METHODS_DIR="$HOME/data/scRFU_radiation_methods_20260916"
SCRFU_PY="$HOME/data/scRFU_regulatory_20260909/environment/bin/python"
NUMBA_CACHE_DIR=/tmp/scrfu-radiation-numba-cache "$SCRFU_PY" \
  -m manuscript.scripts.rp1_14_prepare \
  --raw 'data/Data to be shared_RP P1-14' \
  --rfu data/RFU_out_RP_P1-14 \
  --out "$SCRFU_METHODS_DIR/prepared/rp1_14_v1"
NUMBA_CACHE_DIR=/tmp/scrfu-radiation-numba-cache MPLBACKEND=Agg "$SCRFU_PY" \
  -m manuscript.scripts.rp1_14_longitudinal \
  --prepared "$SCRFU_METHODS_DIR/prepared/rp1_14_v1" \
  --out "$SCRFU_METHODS_DIR/results/rp1_14_v1"
```

Completed output reuse checks scientific input/configuration/code identity and
output hashes. See the active execution state and manuscript for endpoint
results. Legacy eRFU/log10sum scores and mixed-cohort embeddings are not inputs.
