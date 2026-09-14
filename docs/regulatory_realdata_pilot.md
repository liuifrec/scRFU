# RFU-QTL × CD4 regulatory pilot — 2026-09-10 checkpoint

The historical checkpoint below is preserved. **The completed, critically reviewed
2026-09-14 application follows at the end of this document.** Its interpretation
supersedes the provisional decision and incomplete-caQTL statements here.

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

## Completed application and critical review — 2026-09-14

### 1. Dataset integrity and reproducibility

Both pinned archives are now complete: eQTL release
[18261456](https://zenodo.org/records/18261456), 1,743,867,610 bytes, MD5
`e06e21a30576e6e271d974f922793ea6`; caQTL release
[18317808](https://zenodo.org/records/18317808), **4,273,913,755 bytes**, MD5
`0a247f926008e7e7792a7689d327709a`. The original resumable caQTL transfer completed;
separate bounded tail-range downloads were not needed for assembly. No omnibus,
large caQTL log-BF archive or ChromBPNet model was downloaded.

The caQTL inventory has 46 members (including directory entries). Selected chr7
and chr6 nominal members contain 13,191,608 and 20,102,870 rows; conditional files
contain 3,904 and 4,578 rows. All requested chr6/7 nominal and conditional exports
were processed. The previously recovered chr7 nominal and chr6 conditional members
have identical SHA256 values in the verified archive. Full member names, sizes,
hashes and inventories are in `prepared/matos_{eqtl,caqtl}_archive_manifest.json`.

RfuWAS Data 1 remains 2,083 associations, 623 GRCh38 variants, 59 RFUs; all 623
earlier reference checks were reused. Workbook SHA256 remains
`42742f4a30548c1184f8de022a2adc7b6b9abd321e81a01ea740fad5527891d9`.
No allele reversal, inferred effect allele, liftover or cross-study direction
comparison was performed. Correct REF identity does not establish beta orientation.

The persistent Python executable is `environment/bin/python` beneath the external
workspace (Python 3.10). It imports the working checkout's `src/scrfu/__init__.py`.
Versions: pandas 2.3.3, NumPy 2.2.6, pyarrow 21.0.0, openpyxl 3.1.5, scikit-learn
1.7.2; `/usr/bin/Rscript` reports 4.1.2. `logs/environment_freeze.txt` records the
full environment; source-code hashes, input hashes and analysis settings are saved
externally. The environment inherits the existing `rfu_manuscript` environment;
it is persistent, not a fully isolated container. No core dependency was added.

`results/final_pilot/` preserves the complete API outputs and dossiers;
`results/partial_caqtl_checkpoint/` remains unchanged. Completed subprocess exit
statuses and provenance confirmed that the initial postprocessing finished,
including figures and receptor joins. The review added explicit `completion.json`
markers: running before work, complete only after required outputs are hashed.
An early `evidence_counts.json` alone is never evidence of completion.

### 2. Exact association coverage and conditional selection

The actual `scrfu.tl.regulatory_triangulation()` ran with
`allow_allele_reversal=False`; every direction-concordance field remains unknown.

| Evidence | Variants / 623 | RFU associations / 2,083 | RFUs / 59 | Molecular records |
| --- | ---: | ---: | ---: | ---: |
| Represented in eQTL export | 463 | 1,561 | 42 | 39,704 variant–gene tests |
| eQTL query-set Bonferroni support | 446 | 1,544 | 41 | See classified lookup |
| Source conditional eQTL selection | 17 | 120 | 39 | 19 target/rank records, 15 genes |
| Represented in caQTL export | 463 | 1,561 | 42 | 52,920 variant–peak tests |
| caQTL query-set Bonferroni support | 462 | 1,560 | 42 | See classified lookup |
| Source conditional caQTL selection | 32 | 160 | 39 | 41 target/rank records, 35 peaks |
| Both conditional layers at the same variant | 6 | 62 | 36 | 6 eQTL + 10 caQTL target/rank records |

The unchanged eQTL counts reproduce the earlier checkpoint exactly. Complete
caQTL coverage adds 1,167 chr6 tests and seven represented/query-supported variants;
it adds RFU 2811 to the query-supported RFUs. The full predefined query family is
**92,624**, replacing the historical 91,457; its threshold is
**5.398168941095181e-7**. Both-layer query support remains 446 variants / 41 RFUs.

Across either conditional layer there are **43 variants, 218 variant–RFU pairs,
60 within-target selection records and 39 RFUs**. These are not 60 mutually
independent mechanisms. All lie in the TRB neighborhood; no LD calculation was
used to count independent loci. Source signal identifiers include release,
context, file, target and rank, so ranks reused across targets do not merge signals.
Sparse caQTL ranks such as 44, 61 and 80 are preserved. TensorQTL's backward
selection can retain sparse forward indices; rank 80 is not evidence of 80 retained
signals for that peak. The exact publication-run software environment is unavailable.

Thirty-nine RFUs have some conditional eQTL variant and some conditional caQTL
variant. Only 36 have both at the *same* variant. RFUs **1466, 3490 and 4699**
have the two layers only at different variants. Their evidence must not be called
same-variant overlap.

### 3. Source significance, fine mapping and colocalization audit

The pinned author code is
[`4b50d55`](https://github.com/marlmatos/cd4t-qtl-map/tree/4b50d55a9bf8339f3625a4de6da9bf618c288b48).
Nominal mapping writes cis-pair statistics rather than a significance-only table,
subject to study MAF/cis-window/phenotype inclusion. Nevertheless, absence from the
released exports is described as **not represented**, not an established untested
or null association. Neither summary archive supplies the permutation/FDR tables
needed to reconstruct exact marginal significance. `source_marginal_significant`
therefore remains unknown even when conditional selection is observed.

The [supplement workbook](https://zenodo.org/records/18408393) is 67,657,102 bytes;
MD5 `2a8ed13f379946c1a7efce117c961201`, verified independently of its enclosing ZIP.
The cached workbook audit records sheet headers, row counts, titles and scope:

| Table | Content | Relevant exact RFU-variant rows |
| --- | --- | ---: |
| 1 | 159,908 ChromBPNet-scored variants | 39; one author IPS-significant variant |
| 2 | Motif/TF class summary | Not variant-level molecular inference |
| 3 | 5,295 disrupted-motif rows | 2 at `7:142701355:G:C`, ETS1 / IRF4_IRF8 |
| 4 | GWAS study catalog | Not RFU association evidence |
| 5 | 8,423 positive QTL–GWAS colocalization rows | 0 exact RFU-QTL variants |
| 6 | 510 ChromBPNet/GWAS credible-set summaries | 0 exact ChromBPNet RFU variants |

The one IPS-significant variant is a predictive sequence-effect annotation, not
new molecular fine-mapping support. Table 5's `PP.H4.abf` is **QTL–GWAS**, not a
direct gene–peak posterior; `coloc_class=both` does not justify relabeling it.
No standalone molecular gene–peak signal-pair result or candidate molecular CS/PIP
was available in these inspected releases. Molecular colocalization at the RFU
loci is **unresolved**, not disproved by absent positive-table rows.

The author's eQTL fine-mapping script requests coverage 0.95, falls back to 0.10
on failure/no sets, and contains active interactive path/N overrides. Without
per-output run manifests, requested/achieved coverage cannot be assigned to our
candidates. This code snapshot alone does not invalidate published results.
Although an upstream script creates r-squared LD, the consuming eQTL R function
regenerates signed `--r square`, matches alleles using `snp_match`, intersects
variant IDs and orders summary statistics to the LD matrix. Thus it is not
defensible to conclude that r-squared was necessarily passed to SuSiE.

The author's molecular coloc script uses shared aligned variant IDs, at least
200 common variants, `coloc.bf_bf(..., overlap.min=1)` and retains PP.H4 > 0.5.
Its `coloc_results/ca_eqtl_coloc/chr*_coloc_results.csv` and the merged
variant/credible-set tables are the compact outputs to request. Their signal IDs,
coverage and exact RFU-variant membership are needed before reusing any call.
Nearby gene–peak coloc would still not constitute RFU colocalization.

### 4. Candidate mechanisms and the explicitly questioned record

Tier A means the *same RFU-associated variant* was selected conditionally in both
molecular layers. Tier B requires one conditional layer plus query-threshold
support in the other. These are descriptive association tiers, not probabilities
of a mechanism. There are **6 Tier A variants / 62 variant–RFU pairs / 36 RFUs**,
and **36 Tier B variants / 155 pairs / 39 RFUs**; RFU sets overlap.

| Variant | Conditional eGene(s) | Conditional caPeak suffix(es) | Interpretation |
| --- | --- | --- | --- |
| `7:142354123:C:T` | TRBV6-2 | 83324a, 83326a | A; direct receptor architecture |
| `7:142355657:A:G` | TRBV7-2 | 83326a | A; RFU 2905 representative retained |
| `7:142377274:C:CT` | TRBV4-2 | 83315, 83318b | A; receptor-locus indel |
| `7:142713241:G:A` | TRBV28 | 83378a, 83382f, 83383a | A; RFU 1415 retained with revised peaks |
| `7:142812489:C:T` | TRBV30 | 83386bb | A; competing receptor target near EPHB6 region |
| `7:142820419:G:A` | ENSG00000289938 | 83343 | A; non-TCR lncRNA target, no demonstrated gene–peak mechanism |
| `7:142847296:C:T` | EPHB6, rank 1 | None | B; RFU 947 / 83386ba chain downgraded |
| `7:142809879:G:A` | EPHB6 rank 2; ENSG00000289938 rank 2 | None | B; RFU 519, competing targets retained |
| `7:142673687:G:T` | ENSG00000288882 | None | B; second lncRNA association retained |

The four representative RFUs have these *full conditional eGene sets*: 947 and
519 each EPHB6 / ENSG00000289938 / TRBV30; 1415 TRBV28 / TRBV18 / TRBC1 / TRBV7-7 /
ENSG00000288882; 2905 TRBV7-2 / TRBV6-2 / TRBV4-2 / TRBV5-1 / TRBV10-2. All
competing query-supported genes/peaks and source statistics are retained in the
external dossiers. RFU 1415's earlier peak 83383b remains a marginal association
at `142713241`; its conditional lead is another variant. RFU 2905's caQTL for
83326a has source rank 80 and conditional p=4.7431e-20; the original eQTL–caQTL
coincidence is now source-conditional in both layers, not demonstrated coloc.

**`7:142543885:C:T` / peak 83318b is one subset, not the whole caQTL result.**
It occurs in the verified chr7 conditional table as `7:142543885[b38]C,T`, rank 2,
all-CD4 context, conditional beta=-0.0950295, SE=0.023555344,
p=6.868688e-5, permutation p=0.01559844, beta-approximation p=0.01764826.
Data 1 links it only to **RFU 1415** (RFU beta=-0.4750321, p=1.702334e-12).
Effect-allele agreement is unknown. The same peak's *marginal* p is 6.621255e-4;
the strongest marginal caPeak is instead 83361 (p=3.092304e-32).
There are 87 represented eGene tests, none conditional or query-threshold
supported; the minimum is TRBV28 p=4.257914e-5. Consequently this record is the
single conditional-only exception in the residual Tier C category. It is retained
as caQTL evidence but excluded from the Tier A/B shortlist. Tier C does not mean
that this particular source conditional record is false or absent.

### 5. TRB architecture and descriptive locus background

The fixed authoritative interval is [NCBI Gene 6957](https://www.ncbi.nlm.nih.gov/gene/6957),
GRCh38.p14 / NC_000007.14 **142,299,011–142,813,287**, one-based inclusive.
Only 1,059 Data 1 associations are inside that strict interval; the 2,076 chr7
associations occupy the broader TRB neighborhood. EPHB6-associated `142847296`
is flanking, not inside the strict interval.

The conditional eQTL targets comprise **12 direct TCR genes / 14 target-rank
records**, plus **EPHB6 and two lncRNAs / 5 records**. The lncRNAs are annotated
as ENSG00000289938 (142812556–142853827) and ENSG00000288882
(142716831–142719234) in the saved Ensembl GRCh38 lookups. Across 35 conditional
caPeaks, 32 inferred anchors are within TRB, one is within the fixed 250-kb flank,
and two are outside that flank. Anchors are reconstructed from source position
minus `start_distance`; they are **not full peak intervals**. Exact BED intervals
and publication annotation versions remain needed for interval-level interpretation.

RFU identity, V-gene usage, local haplotypes and complex receptor-locus mapping
can all contribute. Structural variation and cross-mapping have not been ruled
out. A non-TCR target label is insufficient to establish a non-TCR mechanism.

The exploratory Matos-only background uses the fixed TRB ±250-kb interval,
100-kb position bins, 0.1-wide MAF bins, SNV/indel class, tested-target-count bins
and minimum target-distance bins. Controls exclude the RFU query set; at least
five controls are required per stratum. It is deterministic standardization,
not a calibrated enrichment test or matched-background simulation.

| Layer | Query/control variants in window | Matched query variants | Query/control threshold-pass rate | Query/control mean −log10(min p) | Query/control conditional rate |
| --- | --- | ---: | --- | --- | --- |
| eQTL | 456 / 1,168 | 76 | 86.84% / 41.87% | 12.10 / 9.60 | 6.58% / 0.032% |
| caQTL | 456 / 1,167 | 94 | 98.94% / 97.98% | 25.87 / 35.93 | 9.57% / 14.20% |

Only 17–21% of represented query variants have this matched support; extrapolating
to the full set is unjustified. caQTL support is saturated in this locus and
controls have stronger minimum-p summaries. Some descriptive eQTL differences
remain, but incomplete RfuWAS ascertainment and LD dependence prevent calibrated
enrichment p-values. The 96–99% raw overlap is not evidence of general enrichment.

### 6. RFU identity beyond TRBV: held-out audit

The upstream reference checkout is `ff7a1ea6aca28444a7ff8ec7ef09745e7341f99d`;
`km5000noMax.Rdata` SHA256 is
`64783074360edeca84b3f49ab9d682263e954752fc09bb208edaabde0fd82553`.
Cached assignment manifests also match RFU.R and trimer-reference hashes.
Assignments used standard `AssignRFUs`, correlation threshold 0.6 and one-based
centroid indices. The reference object itself has no observed V calls. Its exact
checksum compatibility with the RfuWAS reference has not been independently
established; numeric candidate cross-links are therefore conditional on that identity.

Existing [GSE190905](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190905)
provides 27,655 paired receptor records from **six donors**. The final source
paper reports seven donors, so this is explicitly the cached six-donor release,
not the entire published cohort. There are 21,657 eligible threshold-qualified
cells and 12,714 donor V/J/CDR3 amino-acid clonotypes. Source expression clusters
use markers including CCR7/SELL, FOXP3 and cytotoxic genes, independently of RFUs;
exclusion of receptor genes from the original variable-gene list is undocumented.

The frozen model is one-hot TRBV + TRBJ + CDR3 length, with or without RFU,
multinomial logistic regression C=1, no tuning search. Both models use identical
evaluation cells and held-out donors. Encoders, class support and RFU support are
fitted on training donors only. RFUs need 10 training donor-clonotypes and three
donors; rare/unseen RFUs pool into `rare`. Unseen V/J/length categories encode as
zero for that feature. All 13 source classes pass the training support rule in
each actual fold; no evaluation cells were removed by class filtering.

The original analysis additionally excludes training V/J/CDR3 identities present
in the held-out donor. This is **unseen-receptor donor generalization**, stricter
than ordinary donor generalization. A separately labeled sensitivity retains
shared receptors across distinct donors, leaving all other settings unchanged.
The original result is preserved and was reproduced exactly in a probability-
audited rerun: maximum fold-score difference **0.0**. Each audited run saves
86,628 prediction records, explicit class probabilities, evaluation weights,
fold membership hashes and a completion manifest. Recomputing log loss from
saved predictions agrees within 1e-12, with matched cells/labels/weights verified.

**Delta = extended log loss − baseline log loss; positive means worse.**
The same tested helper supplies paired tables and the figure. Primary weights
give each donor clonotype total weight one and balance training donors. Clones
observed in multiple states contribute fractional outcome weights, not duplicate
independent observations. Cell weighting is a clonal-expansion sensitivity.

| Held-out donor | Original clone-weighted baseline | +RFU | Delta | Ordinary donor-generalization delta |
| --- | ---: | ---: | ---: | ---: |
| BFY | 1.766806 | 1.773044 | +0.006238 | +0.006254 |
| RGL | 1.830608 | 1.835367 | +0.004759 | +0.004865 |
| WTS | 1.764356 | 1.769160 | +0.004803 | +0.004841 |
| WZF | 1.942195 | 1.949475 | +0.007280 | +0.007656 |
| XZH | 1.776478 | 1.778449 | +0.001971 | +0.000774 |
| ZBC | 2.035630 | 2.036484 | +0.000854 | +0.000906 |
| Equal-donor mean | **1.852679** | **1.856996** | **+0.004317** | **+0.004216** |

Ordinary donor-generalization mean baseline/+RFU is 1.850996/1.855212. Mean
cell-weighted baseline/+RFU is 2.286472/2.309457 in the original (+0.022985) and
2.244878/2.255544 in ordinary donor generalization (+0.010666). Per-donor scores
for both weighting modes are saved in `paired_donor_deltas.tsv`. No cell-based
confidence interval or calibrated p-value is reported: six donor folds with
overlapping training sets support a descriptive paired assessment, not precise
population-wide uncertainty.

Secondary clone-weighted balanced accuracy changed only slightly: 0.083503 to
0.083727 in the original design and 0.083931 to 0.084436 in ordinary donor
generalization. The predefined primary endpoint, log loss, worsened in both.

Supported RFU categories cover only 0.7–9.4% of original held-out cells. Among
39 conditional-supported RFUs, 28 have any eligible observation and only eight
have observations in at least three donors. RFU 947 has one donor-clone with
TRBV30; 519 and 1415 have none; 2905 has six donor-clones in one donor spanning
TRBV3-1/14/4-1/4-2/4-3. These are inadequate for candidate-specific claims.
In particular the published TRBV12-5 observation for RFU 1415 cannot be tested
here. Matos's queried TRBV12-5 minimum p=0.000292 does not pass the query threshold;
TRBV28 is not an exclusive explanation of that RFU.

**Conclusion: no improvement demonstrated in this dataset/design.** Both
estimands show a small deterioration in clone-weighted log loss across all six
donors. This does not establish that RFUs lack biological value generally, or
that TRBV alone explains every regulatory association. Candidate-specific
validation and independent-cohort replication remain insufficient. Other local
GSE157007 caches supply receptor/sample metadata, not prepared independent
cell-state labels suitable for immediate replication.

### 7. Cell-state and disease annotations

Of 39 Tier A/B RFUs, CD4 annotation is present for 2/39, CD8 for 2/39, TN for
8/39 and Treg for 3/39; 35 lack CD4/CD8 annotation and 28 lack state annotation.
The corresponding all-59-RFU denominators are 2/59 CD4, 13/59 CD8, 9/59 TN and
3/59 Treg. These sparse, nonrandom annotations do not establish CD4/TN/Treg
enrichment. A germline regulatory association measured in CD4 cells does not
contradict a CD8-enriched RFU.

All original 31 disease annotations for 11 RFUs are retained. Complete query-level
coverage adds RFU 2811 with two annotations, yielding 33 for 12 RFUs; high-confidence
Tier A/B retains the original 11 disease-annotated RFUs. Every linked phecode,
effect and p-value is stored in the candidate dossiers, not only celiac disease.
Examples: RFU 947 celiac p=7.0364e-69; RFU 519 celiac p=2.3757e-50; RFU 2905
celiac p=4.3544e-43. These are **downstream genetically predicted RfuWAS phenotype
annotations**, not variant-level GWAS evidence or independent experimental validation.
They do not affect the conditional-evidence tiers. RFU 1415's Treg annotation
does not demonstrate genetic mediation of that state.

### 8. Formal-analysis readiness and author data-request draft

| Required input | Current status |
| --- | --- |
| Dense regional RFU beta/SE or sufficient summary statistics | Missing; Data 1 is significant-only |
| RFU tested-variant universe and missingness | Not released in the inspected repository; scripts reference local matrices/pvar files |
| Verified effect alleles/scales and per-trait sample information | Incomplete |
| Population-compatible allele-aligned LD for each study | Missing for this formal analysis |
| Molecular marginal and conditional statistics | Available at queried variants; full chr6/7 nominal exports retained |
| Source fine-mapping membership, coverage and run provenance | Missing for candidate interpretation |
| Exact RFU reference identity in both studies | Publication-specific checksum unverified |

Draft request, **not sent**: Please provide shareable summary-level RFU-QTL
statistics, including nonsignificant variants, for the GRCh38 TRB neighborhood
and all associated RFUs; the tested variant list, effect/non-effect alleles,
beta/SE and scale, per-trait sample counts, ancestry/covariates, missingness and
RFU normalization/reference checksum. For Matos, please provide candidate-region
gene–peak coloc signal-pair outputs and posterior criteria, credible-set membership
and achieved coverage, PIP, analysis-version/run metadata, peak BED annotations,
and a suitable shareable signed/allele-aligned LD resource or instructions.
We are requesting summary-level data, not participant-level records.

GRCh37 lasso weights are predictive coefficients and cannot fill these gaps.
No new coloc algorithm, formal RFU colocalization, mediation or causal analysis
was run. Source molecular coloc, if supplied, still requires verified RFU variant
membership and cannot be promoted transitively to RFU colocalization.

### 9. Manuscript-ready text

**Results draft.** We applied scRFU regulatory triangulation to 2,083 published
RFU-QTL associations and checksum-verified CD4 molecular-QTL releases accompanying
the Matos preprint. Among 623 GRCh38 RFU-associated variants, 17 overlapped
source-conditional eQTL selections and 32 overlapped conditional caQTL selections;
six overlapped both layers at the same variant, representing 36 RFUs. These
overlaps comprised within-target selection records rather than independent
mechanisms. Most eGene targets were TCR genes, with EPHB6 and two lncRNA targets
providing alternative local regulatory hypotheses. A Matos-only locus comparison
showed nearly saturated caQTL support among matched controls, limiting the
interpretation of broad nominal overlap. In a separate six-donor paired receptor/
expression application, adding RFU categories to TRBV/TRBJ/CDR3-length models
did not improve held-out clone-weighted log loss under either ordinary donor
generalization or exclusion of shared receptor identities. scRFU therefore
organizes traceable regulatory hypotheses here, while added representation–state
information and shared RFU regulatory signals remain unestablished.

**Methods draft.** We preserved published one-based RFU labels, exact GRCh38
chromosome/position/REF/ALT identity and unresolved effect orientation. Verified
chr6/7 molecular statistics were streamed and filtered to the RFU variant set.
A single Bonferroni family comprised 92,624 unique retrieved molecular tests
before RFU/annotation joins. Conditional selections were retained separately with
target-scoped source ranks; unrecoverable marginal FDR and fine-mapping fields
remained unknown. Descriptive locus controls were stratified by position, MAF,
variant class, target opportunity and distance, without calibrated enrichment
p-values. Frozen RFU assignments in GSE190905 were evaluated using identical
leave-one-donor-out folds, training-only preprocessing, fixed logistic models,
donor-balanced training and clone-weighted/cell-weighted evaluation. All held-out
probabilities and paired donor-level scores were saved and audited.

**Limitations draft.** Significant-only RFU statistics, LD dependence, unresolved
effect alleles and fine-mapping provenance, receptor-locus structural complexity,
and missing source molecular-coloc outputs prevent shared-signal or causal claims.
The paired cohort is small, disease/treatment-specific and sparsely covers the
candidate RFUs; the publication-specific RFU reference identity and receptor-gene
influence on source clustering are unresolved. The negative prediction result is
specific to this cohort and fixed design. Disease annotations derive from the
same genetic-prediction framework and are not independent validation.

### 10. Outputs, validation and decision

The generic regulatory module and RFU assignment remain unchanged. Reusable
example extensions implement complete-coverage checks, explicit evidence axes,
dossiers, source-workbook audit, Matos-only descriptive controls and donor-held-out
prediction auditing. New synthetic tests cover missing coverage, workbook integrity,
TRB boundaries, signal strength distinctions, donor/public-clone separation,
clone weighting, paired delta signs, prediction identity/order and incomplete-run
status. Exact final test/lint results and commands are recorded in the execution
state and validation logs.

Final validation: **104 focused tests passed; 466 full-suite tests passed,
4 skipped, 32 warnings**. The final title-only fix additionally passed nine
prediction/completion tests. Ruff lint and format checks passed (218 files
formatted), as did `git diff --check`. Both documentation-format fixes preserve
the Python examples' ASTs. All 21 output hashes in the three completion manifests
were verified after the final figure refresh.

Externally saved figures: `evidence_accounting.pdf`,
`rfu_beyond_trbv_heldout.pdf`, and the earlier complete API coverage figure.
Bulk evidence tables, all phenotype links, model predictions and source files
remain outside Git. The final completion marker covers dossier and figure hashes.

**Decision:** a small set of local association hypotheses merits follow-up,
especially the six same-variant conditional overlaps and competing EPHB6/lncRNA
signals. The current data do not demonstrate RFU colocalization or added cell-state
prediction beyond receptor-gene features. The single most informative next step
is to obtain dense, allele-documented regional RFU-QTL summary statistics with
appropriate study LD, then test shared signals against the competing molecular
targets using an established method. Without those inputs, increasing overlap
counts or selecting attractive disease annotations would add little information.
