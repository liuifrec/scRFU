# RFU-QTL × CD4 regulatory pilot — 2026-09-10 checkpoint

**Decision A, provisionally for the combined analysis:** exact regulatory evidence
justifies locus-level follow-up. The verified eQTL result alone contains 17
independent-QTL variants. The combined result is substantial but remains
**provisional because the caQTL download is incomplete and its archive checksum
cannot yet be verified**. This is not a completed genome-wide regulatory analysis.

All data and generated tables are outside Git at
`~/data/scRFU_regulatory_20260909/`. Existing RFU preparation and reference checks
were reused; the 623-variant reference validation was **not repeated**.

## Inputs and integrity

| Input | Exact file | Bytes | Integrity / scope |
| --- | --- | ---: | --- |
| RfuWAS Supplementary Data 1–4 | `42003_2024_7010_MOESM4_ESM.xlsx` | 831,341 | SHA256 verified during preparation; reused |
| Matos eQTL, release 18261456 | `eQTLs_summary_statistics.tar.gz` | 1,743,867,610 | Full MD5 verified |
| Matos caQTL, release 18317808 | `caQTLs_summary_statistics.tar.gz` | 2,352,865,280 retained of 4,273,913,755 | Incomplete; preserve for resume |
| Matos article supplements | `PMC12870616_supplementary.zip` | 72,278,229 | Read for annotation availability; not used as QTL association input |

Sources: [RfuWAS publication](https://doi.org/10.1038/s42003-024-07010-x),
[Matos eQTL release](https://doi.org/10.5281/zenodo.18261456),
[Matos caQTL release](https://doi.org/10.5281/zenodo.18317808), and
[Matos article supplements](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12870616/supplementaryFiles).

Checksums:

```text
RfuWAS workbook SHA256
42742f4a30548c1184f8de022a2adc7b6b9abd321e81a01ea740fad5527891d9
eQTL archive MD5 (verified)
e06e21a30576e6e271d974f922793ea6
caQTL full archive MD5 (expected, NOT yet verified)
0a247f926008e7e7792a7689d327709a
Article supplement ZIP SHA256 (locally recorded)
f9d7ddb166801b1a930ad55b6f02c039e853545d9d2682fc71a6b813e7115cac
```

The eQTL archive contains 22 nominal Parquet and 22 independent TSV files,
2,461,776,856 uncompressed bytes in total. Only chr6/7 members were extracted:

```text
eQTLs/allcells.cis_qtl_pairs.chr6.parquet
eQTLs/allcells.cis_qtl_pairs.chr7.parquet
eQTLs/allcells.independent_cis_qtl_pairs.chr6.csv
eQTLs/allcells.independent_cis_qtl_pairs.chr7.csv
```

The saved caQTL prefix contains these **complete individual members**:

```text
caQTLs/cd4_qsmooth_cpm_chromatin_narrowpeaks.cis_qtl_pairs.chr7.parquet
  321399762 bytes; SHA256
  534e0e006b8799a62bddf29cf800ef2a6eabd19fc9c798b17190edfe115b9dd2
caQTLs/cd4_qsmooth_cpm_chromatin_narrowpeaks.independent_cis_QTL_pairs.chr6.csv
  1062096 bytes; SHA256
  51ac7b391279d084c8688ed3f5a6490bd0e04f26c8a78cf2ba2a9c526fab6d4f
```

The extracted Parquet is readable and declared member lengths agree. These
locally recorded hashes do **not** replace verification of the publisher's full
archive checksum. Chr6 nominal and chr7 independent caQTL tables remain missing.
Official direct, API, bounded-range, and metadata requests returned 504 after the
initial successful transfers; further identical retries were stopped.

The external workspace occupies approximately **4.6 GiB** at this checkpoint.
No omnibus archive, caQTL log-BF archive, or ChromBPNet models were downloaded.

## RFU-QTL and allele QC

Data 1 supplies **2,083 associations, 623 unique variants, 59 RFUs**, with no
duplicate variant–RFU pairs or unresolved variant identities. All 623/623 REF
sequences, including indels, matched GRCh38 in the previous validation.

There are 616 chr7 variants / 2,076 associations and seven chr6 variants / seven
associations. Of the associations, 1,059 lie within the explicitly defined TRB
interval `chr7:142299011–142813287` (GRCh38, one-based inclusive); all 2,076 chr7
associations fall within its ±1 Mb flanks. Do not describe all of them as lying
inside TRB itself.

Matos identities use `CHR:POS[b38]REF,ALT`. Exact retained identities agree with
the already reference-validated RFU variants; a previous convenience sample also
matched the first Matos allele to REF for 300 eQTL and 300 caQTL variants. This
does not establish the allele used to estimate beta. No effect allele was
invented; `allow_allele_reversal=False`; all output direction-concordance values
are missing. No liftover was performed.

## Evidence definitions

The eQTL lookup streamed **7,668,924 chr6 + 5,007,630 chr7** nominal rows in batches
of 100,000, retaining 39,704 exact variant–gene rows. The recovered chr7 caQTL
lookup streamed **13,191,608** rows, retaining 51,753 variant–peak rows. Repeated
variants across distinct genes/peaks are separate tests, not duplicate variants.

The verified eQTL archive **does not contain the permutation/FDR tables** needed
to recover phenotype-specific significance thresholds, despite their production
in the [Matos mapping code](https://github.com/marlmatos/cd4t-qtl-map).
Therefore source-defined marginal significance is **unknown**, not false.
Independent-QTL membership is separately available and positively identifies
source-selected associations. `rank` is a conditional signal index, not a
credible-set identifier; independent membership does not imply high PIP.

For this exploratory checkpoint, a single analyst-defined Bonferroni family
contains all **91,457 retrieved variant–target tests across both layers**, with
`p < 0.05/91457 = 5.467050089112917e-7`. This was chosen before inspecting their
effect sizes/p-values. It is not the authors' FDR threshold. It controls the
specified retrieved test family, not an unavailable genome-wide RFU-QTL search.
Adding missing chr6 caQTL rows will change the family size and requires reapplying
the filter. Nominal `p < 0.05` is retained only as a descriptive category.

The actual `scrfu.tl.regulatory_triangulation()` input retains rows passing that
Bonferroni threshold **or** matching released independent QTLs. Marginal and
conditional statistics remain separate (`pvalue` versus `conditional_pvalue`,
etc.); independent-only support must not be mistaken for marginal significance.

## Overlap results

Counts are unique RFU-QTL variants, with RFUs in parentheses.

| Category | eQTL chr6/7 | caQTL chr7, provisional |
| --- | ---: | ---: |
| Represented/tested | 463 (42) | 456 (41) |
| At least one nominal p < .05 | 463 (42) | 456 (41) |
| Lookup Bonferroni evidence | 446 (41) | 455 (41) |
| Released independent QTL | **17 (39)** | **Unknown** for chr7 |
| Evidence supplied to API | 446 (41) | 455 (41) |

**446 variants / 41 RFUs have both layers of provisional evidence.** There are
nine caQTL-only variants, zero eQTL-only variants, and 168 variants with no
qualifying evidence in the supplied tables. The last category includes untested
and unavailable data and is not evidence of biological absence.

The verified independent-eQTL result comprises **19 variant–gene pairs, 15 genes,
17 variants, 39 RFUs, and 153 variant–RFU–gene chains**. The recovered chr6
independent caQTL file has zero exact matches to the seven chr6 RFU variants;
this does not establish the total independent-caQTL overlap.

The combined API output contains 4,138 RFU–eQTL matches and 7,367 RFU–caQTL matches,
covering 16 genes and 26 peak IDs. These counts are not independent loci.
Proximity analysis was not triggered: exact support is not sparse.

## Representative evidence chains

Peak suffixes below abbreviate the original `cd4_atac_summits_peak_` prefix.
eQTL p-values in this table are from the **released independent** tables;
caQTL p-values are **marginal and provisional**.

| GRCh38 variant → RFU | Independent gene (p) | Exact caQTL peak (p) | RFU context / strongest disease |
| --- | --- | --- | --- |
| `7:142847296:C:T` → 947 | EPHB6 (3.45e-31) | 83386ba (1.67e-29) | Cell annotation absent; celiac disease, 7.04e-69 |
| `7:142809879:G:A` → 519 | EPHB6, rank 2 (3.77e-21) | 83386bb (3.85e-50) | TN; celiac disease, 2.38e-50 |
| `7:142713241:G:A` → 1415 | TRBV28 (9.65e-77) | 83383b (1.98e-74) | Treg; no Data 4 disease association |
| `7:142355657:A:G` → 2905 | TRBV7-2 (2.96e-28) | 83326a (1.86e-95) | CD8-enriched; celiac disease, 4.35e-43 |

These are overlapping associations, not demonstrated gene–peak regulatory
connections or causal paths. Other genes/peaks can share the same variant.
EPHB6 is an interesting non-TRBV candidate, but its mechanism is not established.
TRBV/TRBC signals may partly reflect inherited receptor-gene usage rather than a
general downstream CD4-state program; LD, structural complexity and read-mapping
effects at the receptor locus need explicit follow-up.

Among the 41 supported RFUs, source Data 2 labels two CD4-enriched (2449, 3766),
two CD8-enriched (2887, 2905), and leaves 37 unannotated. Data 3 labels eight TN
and three Treg, with 30 unannotated. Continuous source metrics are retained;
missing labels are not interpreted as a distinct cell type. CD8 enrichment is
not contradicted merely by observing a shared germline regulatory association
in purified CD4 cells.

Data 4 contributes **31 links across 11 RFUs**. The strongest association is
RFU 947–celiac disease (p=7.0364e-69, source effect −0.0135783); RFU 519–celiac
disease is p=2.3757e-50. RFU 1986 links to rheumatoid arthritis (p=4.0733e-10),
and RFU 267 to disorders of iron metabolism (p=1.8986e-8). These are RFU–phenotype
associations, joined by unchanged published RFU labels, not variant-level GWAS
evidence. Their effect estimates are not compared with molecular-QTL directions.

## Outputs and reproducibility

`results/eqtl_independent_checkpoint/` contains the verified eQTL-only API run.
`results/partial_caqtl_checkpoint/` contains the provisional combined API run:

```text
regulatory_hits.tsv                 rfu_regulatory_summary.tsv
variant_regulatory_summary.tsv      unmatched_variants.tsv
harmonized_rfu_qtl.tsv               provenance.json
matos_eqtl_lookup_classified.tsv     matos_caqtl_lookup_classified.tsv
matos_eqtl_independent_exact.tsv     matos_caqtl_independent_exact.tsv
regulatory_hits_with_cellstate.tsv   regulatory_rfu_disease_links.tsv
ranked_regulatory_candidates.tsv    overlap_counts.json
regulatory_evidence_coverage.png     regulatory_evidence_coverage.pdf
```

The caQTL-independent exact output is chr6-only at this checkpoint. Unknown total
counts are JSON null. Ranking counts observed evidence layers and annotation
availability; it is a lower bound while caQTL independence is missing and is not
a causal probability. Correlated variants and uneven annotation can dominate it.

To reuse the prepared files (no download, reference-QC rerun, or nominal rescan):

```bash
cd ~/Github/scRFU
export PYTHONPATH=/tmp/scrfu-regulatory-deps:$PWD
python examples/run_regulatory_realdata_pilot.py \
  --root ~/data/scRFU_regulatory_20260909 --partial-caqtl
```

`examples/recover_matos_caqtl_checkpoint.py --root ...` reproduces the recovery
only if its outputs need rebuilding. `examples/rfuwas_regulatory_context.py`
reproduces the earlier preparation/reference QC and does not need running now.
The isolated `/tmp/scrfu-regulatory-deps` supplies pyarrow 21.0.0 and openpyxl
3.1.5 (plus et-xmlfile 2.0.0); these are example dependencies, not new scRFU core
dependencies. The analysis itself uses existing pandas/numpy/matplotlib.

Next completion commands, after Zenodo resumes serving the file:

```bash
curl -fL -C - --connect-timeout 20 --speed-limit 1024 --speed-time 120 \
  'https://zenodo.org/records/18317808/files/caQTLs_summary_statistics.tar.gz?download=1' \
  -o ~/data/scRFU_regulatory_20260909/sources/matos/caQTLs_summary_statistics.tar.gz
# If interrupted, rerun this command later: -C - preserves/resumes the prefix.
# Avoid short --max-time plus automatic curl retries that can restart transfers.
python examples/matos_regulatory_lookup.py \
  --root ~/data/scRFU_regulatory_20260909 --layer caqtl
# This verifies the expected full MD5 before extracting/scanning chr6/7.
python examples/run_regulatory_realdata_pilot.py \
  --root ~/data/scRFU_regulatory_20260909
```

The eQTL lookup is complete and need not be rerun. The full run uses the existing
eQTL preparation, adds the missing caQTL rows/independent signals, recalculates the
joint test family, and writes primary results outside the provisional directory.

Validation: **93 focused tests passed; 455 full-suite tests passed, four skipped**
(32 existing dependency warnings). Ruff check, formatting checks for the seven
new Python files, and `git diff --check` passed. Synthetic full/partial pilot runs
exercise the actual API, provenance, annotation joins, and unknown direction /
independence handling. The generic regulatory module was not changed.

## Scientific limits and next step

Data 1 is a significant-only RFU-QTL table, with no genome-wide RFU null/background
statistics. Neither the large overlap count nor the Bonferroni filter establishes
enrichment versus a genomic background. The GRCh37 lasso weights were not used.
There is no LD-aware matching, RFU credible set, formal colocalization, mediation,
causal inference, or antigen-specificity inference here. Independent-QTL status
is not fine-mapping probability. Missing source FDR thresholds, peak interval
annotations and sample-level genotype metadata limit interpretation.

Finish caQTL integrity/coverage first. Then prioritize TRB/TRBV28 and
EPHB6/TRBV30-region loci rather than treating hundreds of correlated variants as
independent discoveries. Formal colocalization needs full locus-level RFU-QTL
and molecular-QTL statistics (not just significant RFU hits), verified effect
alleles and scales, per-study sample sizes/ancestry, compatible variant coverage,
and suitable LD for models allowing multiple signals. Obtain those data and
resolve receptor-locus mapping/structural variation before fitting a separate
validated coloc/SuSiE method. Data 4 disease associations alone cannot supply the
variant-level disease statistics needed for a three-trait analysis.
