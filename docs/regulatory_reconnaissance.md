# Bounded public-resource compatibility audit — 2026-09-08

Audited implementation: `a8526bb`, branch `feature/regulatory-triangulation`.
No package changes or genome-scale joins were made. **A first exact-overlap run is
feasible using genuine public RFU-QTL associations and Matos molecular QTLs in
GRCh38. Direction concordance must remain disabled until allele coding is resolved.**

Sources inspected:

- [Matos eQTL concept record](https://doi.org/10.5281/zenodo.18261455), resolving
  to version **18261456**, dated 2026-01-28.
- [Matos caQTL concept record](https://doi.org/10.5281/zenodo.18317807), resolving
  to version **18317808**, dated 2026-01-28.
- [Matos code](https://github.com/marlmatos/cd4t-qtl-map/tree/4b50d55a9bf8339f3625a4de6da9bf618c288b48).
- [RfuWAS code](https://github.com/YuhaoTan2/RfuWAS/tree/f71578971a9fce671413072d1f1ade09c911c670),
  [paper](https://www.nature.com/articles/s42003-024-07010-x), and actual supplementary workbook.

## Files and transfer bounds

| Version record | Exact downloadable file | Compressed bytes | MD5 from Zenodo |
|---|---|---:|---|
| 18261456 | `eQTLs_summary_statistics.tar.gz` | 1,743,867,610 | `e06e21a30576e6e271d974f922793ea6` |
| 18261456 | `eQTL_lbfs.tar.gz` | 1,381,008,182 | `352b04dffdb3066528b44cabefad350e` |
| 18317808 | `caQTLs_summary_statistics.tar.gz` | 4,273,913,755 | `0a247f926008e7e7792a7689d327709a` |
| 18317808 | `caQTL_lbfs.tar.gz` | 15,307,898,585 | `6ec6a0eb07efee342c955ea9cc145181` |

Only the first **131,072 compressed bytes of each** were transferred via HTTP
Range. These are incomplete archive samples, not verified complete downloads.
Full uncompressed archive sizes and complete member inventories remain unknown;
obtaining them from ordinary tar.gz requires further streaming/decompression.
The following exact members/sizes were read from real tar headers:

| Archive | Observed member | Uncompressed bytes |
|---|---|---:|
| eQTL summary | `eQTLs/allcells.independent_cis_qtl_pairs.chr16.csv` | 76,664 |
| eQTL summary | `eQTLs/allcells.cis_qtl_pairs.chr20.parquet` | 56,111,770 |
| caQTL summary | `caQTLs/cd4_qsmooth_cpm_chromatin_narrowpeaks.independent_cis_QTL_pairs.chr10.csv` | 776,404 |
| eQTL LBF | `eQTL_lbfs/All_CD4T_cells_chr14_lbfs.txt` | 144,400,828 |
| caQTL LBF | `CD4T_chromatin/CD4T_chromatin_chr4_lbfs.txt` | 1,321,580,601 |

The combined 7.8 GB archive was **not downloaded**. ChromBPNet metadata lists
`ChromBPNet Models.zip`, **599,091,007 bytes**; it was **not downloaded** and is not
an association-statistics substitute.

Small complete downloads, retained only under `/tmp`:

- [Supplementary workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs42003-024-07010-x/MediaObjects/42003_2024_7010_MOESM4_ESM.xlsx):
  **831,341 bytes**; SHA256 `42742f4a30548c1184f8de022a2adc7b6b9abd321e81a01ea740fad5527891d9`.
- [Lasso ZIP](https://raw.githubusercontent.com/YuhaoTan2/RfuWAS/f71578971a9fce671413072d1f1ade09c911c670/models/lasso_weights_tsv.zip):
  **1,202,743 bytes**, **1,933,490 uncompressed bytes**; 1,351 TSVs plus directory;
  SHA256 `0b12d0c35cee75903580014b396fb0183cf9168c163ed44e72f3d022ef7c7376`.

Thus biological archive/workbook transfers totaled 2,558,372 bytes, plus small
metadata, source-code and documentation reads. No biological data are committed.

## Actual Matos schemas

The sampled independent tables are **TSV despite their `.csv` suffix**, with a
leading unnamed pandas row-index column. Exact named columns:

```text
phenotype_id num_var beta_shape1 beta_shape2 true_df pval_true_df
variant_id start_distance end_distance ma_samples ma_count af
pval_nominal slope slope_se pval_perm pval_beta rank
```

`phenotype_id` contains gene symbols/Ensembl IDs for eQTL, and chromatin feature
IDs such as `cd4_atac_summits_peak_9601` for caQTL, not interval coordinates.
`rank` identifies conditional signals; it is **not a credible-set ID**.
`slope`, `slope_se`, and `pval_nominal` map to beta, SE, and nominal p-value.
Preserve the permutation/beta-approximation p-values separately.

The [eQTL mapping script](https://github.com/marlmatos/cd4t-qtl-map/blob/4b50d55a9bf8339f3625a4de6da9bf618c288b48/qtl_mapping/eqtl_mapping/008_ciseQTLTensor_allcells.py)
and [caQTL mapping script](https://github.com/marlmatos/cd4t-qtl-map/blob/4b50d55a9bf8339f3625a4de6da9bf618c288b48/qtl_mapping/caqtl_mapping/008_caQTL_narrowpeaks_tensor-1mb.py)
verify tensorQTL nominal, permutation, and conditional-independent workflows.
They write these basename patterns; existence of every chromosome/member has
**not** been established from the bounded archive samples:

```text
eQTLs/allcells.cis_qtl_pairs.chr{N}.parquet
allcells.cis_qtl_pairs.chr{N}.csv
allcells.sig_cis_qtl_pairs.chr{N}.csv
allcells.independent_cis_qtl_pairs.chr{N}.csv
cd4_qsmooth_cpm_chromatin_narrowpeaks.cis_qtl_pairs.chr{N}.parquet
cd4_qsmooth_cpm_chromatin_narrowpeaks.cis_qtl_pairs.chr{N}.csv
cd4_qsmooth_cpm_chromatin_narrowpeaks.sig_cis_QTL_pairs.chr{N}.csv
cd4_qsmooth_cpm_chromatin_narrowpeaks.independent_cis_QTL_pairs.chr{N}.csv
```

The `sig` filename alone is not a significance guarantee: mapping code writes the
entire `map_cis` result after adding q-values; filter those q-values explicitly.
Conditional mapping is gated at FDR 0.05. A first lead-variant run and a later
all-nominal-variant run answer different overlap questions.

The sampled eQTL LBF table has exact columns:

```text
lbf_1 ... lbf_10 variant_id.x PIP region gene nsnps variant_id.ss
variant_id.y start_distance af ma_samples ma_count pval_nominal
beta slope_se pos chr a0 a1
```

The sampled caQTL LBF table replaces `gene` with `peak` and additionally contains
`V1` immediately after `variant_id.ss`. Both are TSV and **already contain PIP**;
no reconstruction from log-BFs is needed. `nsnps` is the number of variants in the
fine-mapping region, not sample size. The ten `lbf_*` columns describe SuSiE
components, not credible-set IDs.

Neither sampled LBF header contains `cs` or `coverage`. Source fine-mapping code
writes separate `*_credible_sets.txt` products with `cs`, `pip`, `region`,
`coverage`, target, and original variant ID, but their inclusion in the released
archives was **not verified**. Do not invent set IDs from LBF column numbers or
from the region. PIP-only support can proceed; credible-set overlap cannot yet.

Context: eQTL all-CD4 pseudobulk, caQTL CD4 ATAC narrow peaks with qsmooth/CPM
processing. Mapping paths identify Ashkenazi genotype subsets and filenames with
363 (scRNA) / 362 (ATAC), but these names do not establish exact analysis sample
counts or per-variant missingness. Full sample/covariate manifests were not
identified in the metadata/samples inspected. Verify model sample sizes before
any formal statistical analysis; no N was inferred from `af`, `ma_count` or `nsnps`.

### Matos build and allele coding

Actual IDs include `16:38552[b38]C,T`, `14:20038319[b38]G,C`, and indels such as
`14:20041347[b38]GA,G`. Numeric chromosome names and one-based GRCh38 positions
are consistent with the core. The ID-generation command uses `REF,ALT`, while
fine-mapping parsing assigns `a0` to the second allele and `a1` to the first and
comments label these reference/effect respectively.

A one-base [UCSC hg38 sequence query](https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr14;start=20038318;end=20038319)
returned **G** for `14:20038319`, supporting first-allele REF for that observed
record. Its released LBF row has `a0=C`, `a1=G`. Thus **do not map a0→REF or
blindly map a1→effect_allele**. One reference check does not validate all variants.

Mapping uses tensorQTL's PLINK reader; its current source flips loaded dosage
(`2 - bed`). Exact historical loader versions and released BIM/counted alleles
remain unverified. A safe first run uses ID `ref-alt` after reference validation,
leaves effect_allele missing, and disables reversed-key matching. For signs,
verify original genotype dosage coding and which beta is rejoined in LBF output;
LD-aligned `variant_id.x` and original `variant_id.ss` must not be interchanged.

## RfuWAS: distinct biological objects

All workbook sheets have a title row, then column headers: pandas `header=1`.

| Artifact | Observed content | Appropriate role |
|---|---|---|
| Workbook `Data 1` | **2,083** rows; `SNP RFU beta t.stat p.value` | **Actual significant genome-wide RFU-QTL associations** |
| Workbook `Data 2` | **5,000** RFUs; CD4/CD8 annotations | RFU-level annotations |
| Workbook `Data 3` | **168** rows; CM/TN/Treg/Tscm annotations | RFU-level state annotations, not all predictable RFUs |
| Workbook `Data 4` | **2,309** rows; predicted RFU–phecode associations | RFU-level phenotype evidence, not variant GWAS |
| `lasso.ukbflip_tsv/RFU*_weights.tsv` | **1,351** fitted RFU model files | Prediction weights / publicly available model RFU list |
| Supplementary Information PDF, Table 2 | Paper identifies significant HLA-haplotype–RFU associations | HLA association evidence; not biallelic variant IDs |

`Data 1` is the correct canonical `rfu_qtl` input. Example:
`7_142581165_C_G`, RFU `1415`, beta `0.81752338744202002`,
t `19.414814941320898`, p `1.04750095161761E-66`.
The [association analysis](https://github.com/YuhaoTan2/RfuWAS/blob/f71578971a9fce671413072d1f1ade09c911c670/scripts/rfuQTL/rfu_eqtl_analyze.R)
explicitly declares **hg38**, exports significant rows, and increments internal
zero-based RFU indices when producing the publication table. **Do not add one
again** to workbook RFU labels. Confirm scRFU's frozen-reference label namespace
before joining existing RFU summaries; numeric `1415` is not automatically the
same string as `RFU1415`.

The [genotype pipeline](https://github.com/YuhaoTan2/RfuWAS/blob/f71578971a9fce671413072d1f1ade09c911c670/scripts/rfuQTL/rfuQTL.sh)
constructs `CHR_POS_REF_ALT` IDs and exports PLINK2 `Av` dosage for MatrixEQTL.
[PLINK2 documentation](https://www.cog-genomics.org/plink/2.0/data#export) says REF
is counted by default. This supports a **REF-counted effect inference from code**,
but the publication workbook omits counted-allele/PLINK-version metadata; retain
unknown effect allele in the first run until that inference is verified against
the original `.traw` or authors. SE is absent; it could be derived as
`abs(beta / t.stat)` under the reported linear-model test, but is unnecessary for
our descriptive overlap and was not manufactured in the smoke test. No RFU-QTL
PIP/credible sets or full genome-wide nonsignificant statistics were identified.

Exact annotation headers:

```text
Data 2: RFU statistic.ms p.value.ms median_ratio.ms enrichment.ms
        statistic.zumla p.value.zumla median_ratio.zumla enrichment.zumla enrichment
Data 3: RFU Friedman_chi2 Friedman_p Enriched_Group CM_median TN_median
        Treg_median Tscm_median adjusted_Friedman_p
Data 4: RFU effect se zscore pvalue n_samples phecode description group
```

Exact lasso header: `#CHROM POS_hg19 RSID REF ALT weight allele_freq`.
For example `lasso.ukbflip_tsv/RFU3340_weights.tsv` has 13 predictors and is
827 uncompressed bytes. Some RSID cells contain comma-separated aliases.
The [export code](https://github.com/YuhaoTan2/RfuWAS/blob/f71578971a9fce671413072d1f1ade09c911c670/scripts/lasso/fit_final_lasso.R)
confirms penalized coefficients and `allele_freq` from `ALT_FREQS`. There are no
association p-values or SEs. The coordinates are explicitly **GRCh37/hg19**.
Weights must never be relabeled as RFU-QTL beta, even after valid liftover.

The HLA PDF is publicly linked as `42003_2024_7010_MOESM2_ESM.pdf`; its table was
not parsed in this bounded audit. The GitHub repository has analysis scripts but
no committed `tables/rfuqtl.tsv` or HLA result tables; the workbook supplies the
former result publicly. Full prediction-validation metrics / complete predictable
RFU annotations were not identified as a separate public table.

## Compatibility matrix

| Resource | File | Build | Variant representation | Effect allele | Target | Direct? / preprocessing |
|---|---|---|---|---|---|---|
| Matos eQTL independent | `allcells.independent_cis_qtl_pairs.chr16.csv` sampled | 38 | `CHR:POS[b38]REF,ALT` | Unresolved | Gene | Existing tensorqtl adapter; TSV; explicit allele order |
| Matos caQTL independent | `…independent_cis_QTL_pairs.chr10.csv` sampled | 38 | Same | Unresolved | Peak ID | Existing tensorqtl adapter; peak annotation needed only for intervals |
| Matos nominal | `allcells.cis_qtl_pairs.chr20.parquet` observed | 38 expected from workflow | Not decoded in this audit | Unresolved | Gene | Parquet→DataFrame/TSV with separate analysis reader; inspect actual columns |
| Matos fine-mapping | `All_CD4T_cells_chr14_lbfs.txt`, `CD4T_chromatin_chr4_lbfs.txt` sampled | 38 | Original `.ss`, LD-aligned `.x`, rejoined `.y` | Unresolved | Gene / peak | Existing susie adapter handles `.ss`, PIP; no CS membership available in sampled headers |
| RfuWAS RFU-QTL | XLSX `Data 1` | **38** | `CHR_POS_REF_ALT` | REF inferred from pipeline; verify | RFU | Read header row 2; rename columns; split ID; preserve published RFU label |
| RfuWAS lasso | `RFU*_weights.tsv` | **37** | Coordinates + REF/ALT + RSID | Model-specific; not audited fully | RFU | **Not RFU-QTL evidence**; separate prediction layer |
| RfuWAS disease | XLSX `Data 4` | Not applicable | No variant | Not applicable | RFU–phecode | Join by RFU; **not** generic variant `gwas` input |

## Audit of a8526bb and smoke result

- No core or adapter bug requiring a patch was demonstrated by these samples.
  The existing explicit allele-order requirement and missing-direction behavior
  are necessary; replacing them with the fine-mapping comments would be unsafe.
- The earlier example's generic filenames were illustrative, not verified archive
  names. Real summary archives also contain **Parquet**, which manifest
  `pd.read_csv` cannot ingest directly; convert outside the core first.
- `.csv` can mean TSV here. Default TSV reading is correct for sampled independent
  files; setting `sep=','` based on the extension would break them.
- Existing `susie` adapter correctly prefers `.ss`, maps PIP without reconstruction,
  and preserves per-target rows. Do not deduplicate solely on variant ID: variants
  recur across genes/peaks. Retain source/context/region and association type.
- Earlier example credible-set filenames and the combined-release label do not
  demonstrate that CS membership is released in the separate archives. Use the
  pinned separate version IDs, and omit credible-set analyses until verified.
- Core accepts opaque RFU IDs with coordinates; it cannot distinguish lasso
  weights from association betas semantically. Correct upstream artifact selection
  is the caller's responsibility. Likewise Data 4 must not become variant GWAS.
- **50 real rows normalized successfully**: ten each from four sampled Matos TSVs
  through `adapt_matos_table`, and ten from RFU-QTL `Data 1` through
  `normalize_regulatory_variants`. All effect alleles were intentionally missing;
  both LBF subsets retained all ten PIPs. This verifies structure, not whole-file
  allele validity, scientific direction, overlap yield, or large-file performance.

Temporary audit material: `/tmp/regulatory-real-smoke.json`,
`/tmp/smoke_regulatory_real.py`, `/tmp/inspect_rfuwas_xlsx.py`, workbook, ZIP and
four range samples. `/tmp` may be cleared; the commands below identify the sources.
The workbook was inspected using Python ZIP/XML because openpyxl is absent from
the current environment. No packages were installed and no tests/code were changed.

## First run after reset

Start with RFU-QTL Data 1 × eQTL/caQTL **conditional lead variants**, using GRCh38
on both sides; no liftover is required for that choice. This is a limited
lead-variant overlap analysis, not full locus-level overlap or colocalization.
First inventory the downloaded archives and confirm chr6/chr7 member names; do not
assume they equal the script patterns. These commands are for the next session,
**not executed in this audit**:

```bash
mkdir -p /tmp/scrfu-first-real
curl -fL --retry 2 'https://zenodo.org/api/records/18261456/files/eQTLs_summary_statistics.tar.gz/content' -o /tmp/scrfu-first-real/eQTLs_summary_statistics.tar.gz
curl -fL --retry 2 'https://zenodo.org/api/records/18317808/files/caQTLs_summary_statistics.tar.gz/content' -o /tmp/scrfu-first-real/caQTLs_summary_statistics.tar.gz
curl -fL 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs42003-024-07010-x/MediaObjects/42003_2024_7010_MOESM4_ESM.xlsx' -o /tmp/scrfu-first-real/rfuwas.xlsx
md5sum /tmp/scrfu-first-real/*tar.gz
tar -tzf /tmp/scrfu-first-real/eQTLs_summary_statistics.tar.gz > /tmp/scrfu-first-real/eqtl_members.txt
tar -tzf /tmp/scrfu-first-real/caQTLs_summary_statistics.tar.gz > /tmp/scrfu-first-real/caqtl_members.txt
rg 'independent.*chr(6|7)\.' /tmp/scrfu-first-real/*members.txt
```

These two summary archives require **6.02 GB compressed** in total; confirm disk
and transfer budget then. Extract only confirmed chr6/chr7 independent members
with `tar -xzf ARCHIVE -C OUTPUT EXACT_MEMBER`. Retain tar checksums and source
member paths in provenance. Do not download the 15.3 GB caQTL LBF archive for the
first run. For later eQTL PIP support, the exact URL is:
`https://zenodo.org/api/records/18261456/files/eQTL_lbfs.tar.gz/content`.
Ask for selective chromosome exports/indexed access before taking on caQTL LBFs.

In a separate analysis environment with openpyxl available (or using the existing
standard-library workbook extraction), prepare Data 1 as follows:

```python
rfu = pd.read_excel('/tmp/scrfu-first-real/rfuwas.xlsx', sheet_name='Data 1', header=1)
rfu = rfu.rename(columns={'RFU': 'rfu_label', 'p.value': 'pvalue'})
rfu['variant_id'] = rfu['SNP'].str.replace('_', ':', regex=False)
rfu['genome_build'] = 'GRCh38'
rfu['source'] = 'RfuWAS Supplementary Data 1'
rfu['release'] = 's42003-024-07010-x'
# Leave effect_allele missing until verified. Published RFU labels are already one-based.
```

Use `adapt_matos_table(..., format='tensorqtl', allele_order='ref-alt',
genome_build='GRCh38', source='Matos-CD4', release='18261456' or '18317808')`
after validating REF against GRCh38. Preserve conditional `rank` and explicitly
label the input analysis as conditional-independent. Then call the core directly:

```python
result = scrfu.tl.regulatory_triangulation(
    rfu, eqtl=eqtl, caqtl=caqtl, allow_allele_reversal=False,
    input_metadata=metadata_with_checksums_contexts_and_filters,
)
```

The current manifest CLI uses the core's reversal-enabled default, so the direct
API is preferable for this stricter first run. No positive overlap is promised;
report unmatched variants and ascertainment differences. Further prerequisites:
verify historical effect-allele coding for sign comparisons, recover actual CS
membership if needed, and verify RFU reference/label identity before adding
phenotype or single-cell evidence. Public significant-only RFU-QTL data are
insufficient for a formal colocalization analysis.
