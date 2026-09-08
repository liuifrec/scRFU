# RFU regulatory triangulation and genetic evidence

This experimental, offline layer connects **supplied** RFU-QTL associations to
molecular QTL and GWAS evidence. It does not estimate QTLs or change RFU assignment.
A useful working hypothesis is genotype → regulatory mechanism → T-cell state →
RFU → phenotype/disease. These arrows describe questions for follow-up, not
relationships established by this analysis.

For example: “RFU-QTL variant X has overlapping CD4 caQTL, eQTL, and GWAS evidence,
nominating gene Z as a candidate regulatory link.” Exact variant overlap alone
does not establish that the gene regulates the peak, that either mediates the
RFU association, or that the RFU mediates disease risk.

## Supported evidence and canonical schemas

All inputs are pandas DataFrames. `scrfu.tl.regulatory_evidence_schema(layer)`
returns a `RegulatoryEvidenceSchema` with required/optional columns for
`rfu_qtl`, `eqtl`, `caqtl`, or `gwas`.

| Layer | Required target | Interpretation |
|---|---|---|
| `rfu_qtl` | `rfu_label` | Supplied variant–RFU abundance association |
| `eqtl` | `gene` | Supplied variant–gene expression association |
| `caqtl` | `peak` | Supplied variant–chromatin feature association |
| `gwas` | `trait` | Supplied variant–trait association |

Each nonempty row requires `source` and `genome_build`. Normally it also supplies
`chromosome`, `position`, `ref`, and `alt`. Alternatively, a `variant_id` of the
form `1:12345:A:G` supplies these coordinates. rsID-only or incomplete-coordinate
rows are retained with unresolved identity and reported in QC; they cannot match.
A schema-less empty DataFrame is accepted. Missing target/source/build values,
invalid supplied statistics, or conflicting coordinate IDs raise `ValueError`.

| Common columns | Contract |
|---|---|
| `variant_id`, `rsid` | Original identifier and optional `rs123` alias; rsIDs are never sufficient matching keys |
| `chromosome`, `position`, `ref`, `alt` | One-based genomic position; forward-genomic allele sequences |
| `genome_build` | Required assembly, canonicalized for hg19/GRCh37 and hg38/GRCh38 |
| `beta` | Optional signed effect on the source outcome scale |
| `se` | Optional positive standard error on that effect scale |
| `pvalue` | Optional probability in [0, 1]; unadjusted/adjusted semantics belong in source metadata |
| `allele_frequency` | Optional frequency in [0, 1]; record its allele convention in metadata |
| `odds_ratio` | Optional positive GWAS OR; do not supply both beta and OR in one row |
| `effect_allele` | Optional explicit allele to which beta/OR refers; must equal REF or ALT |
| `strand` | Optional `+`/`-` orientation declaration; negative strand direction is unsupported |
| `pip` | Optional supplied posterior inclusion probability in [0, 1] |
| `credible_set_id`, `locus_id` | Optional set and locus identifiers; locus IDs must disambiguate reused set IDs |
| `context` | Optional cell type, assay context, condition, or analysis stratum |
| `source`, `release`, `url`, `sha256` | Source name, version, resource URL, optional file checksum |

`beta`, `se`, and `pvalue` may be absent, including in credible-set-only exports.
Canonical columns have exact names: explicitly rename source `effect_size` to
`beta`, for example, after checking its definition. Unknown columns are retained;
this permits coverage, sample size, ancestry, study phenotype coding, assay,
RFU-reference ID, and upstream QC annotations. Scalars are expected in cells.
Peak strings may be intervals or feature IDs; the core does not parse intervals,
map genes to peaks, or perform interval overlaps.

Normalize directly, or provide per-layer defaults through `input_metadata`:

```python
import pandas as pd
import scrfu

rfu_qtl = pd.DataFrame({
    "variant_id": ["chr1:12345:a:g"], "genome_build": ["hg38"],
    "rfu_label": ["RFU1"], "beta": [0.3], "effect_allele": ["G"],
    "source": ["my-rfuQTL"], "release": ["analysis-v1"],
})
eqtl = pd.DataFrame({
    "variant_id": ["1:12345:A:G"], "genome_build": ["GRCh38"],
    "gene": ["GENE1"], "beta": [0.2], "effect_allele": ["G"],
    "pip": [0.98], "context": ["CD4"],
    "source": ["local-eQTL"], "release": ["release-v1"],
})
result = scrfu.tl.regulatory_triangulation(rfu_qtl, eqtl=eqtl)
print(result.rfu_summary)
```

This illustration is synthetic. No biological annotation or data are bundled.

## Variant and effect harmonization

`normalize_regulatory_variants` copies its input. It strips chromosome `chr`
prefixes, normalizes numeric chromosomes and X/Y/MT, accepts positive integer
positions (including `"12345.0"`), and uppercases A/C/G/T allele sequences.
`variant_key` is `CHROM:POS:REF:ALT`; **use it together with `genome_build`**.
Only hg19/GRCh37 and hg38/GRCh38 are assembly aliases. Other assembly names are
retained exactly; different names require external verification. Missing or mixed
builds fail before matching, even for variants on disjoint chromosomes. Patch
assembly names are not silently mapped to their parent assembly.

The core does not validate REF against a FASTA, query dbSNP, split multiallelic
records, complement alleles, normalize/left-align indels, expand LD proxies, or
perform liftover. Supply externally normalized, biallelic coordinates. Anchored
indels consisting of A/C/G/T sequences are accepted, but equivalent differently
represented indels will not match. Unresolved rows remain in RFU summaries and QC,
and are excluded from per-variant summaries and credible-set intersections.

By default `allow_allele_reversal=True` additionally matches A:G to G:A at the
same position and build. Such rows have `alleles_reversed=True` and
`shared_variant_exact=False`. Supplied REF/ALT values and original effects remain
in output, making orientation discrepancies reviewable. Set this option to
`False` to require identical ordered allele keys throughout triangulation.

Direction comparison requires explicit `effect_allele` on both records. ALT is
never inferred as the effect allele. The external beta is sign-flipped only when
its effect allele differs from the RFU effect allele. GWAS OR is converted to
log(OR) for this sign comparison. A/T and C/G pairs require explicit `strand="+"`
on both records, even for an exact key; negative strand direction is unsupported.
Missing alleles, unresolved strand, absent effects, and zero effects leave
`direction_concordant` missing with an explanatory `direction_status`.
`allele_harmonized=True` means alleles could be oriented; effects may still be
missing. Concordance only describes effect **signs** on separately defined outcome
scales; it does not compare effect magnitudes or demonstrate mediation.

## Result tables and evidence tiers

`RegulatoryTriangulationResult` contains:

| Attribute | Contents |
|---|---|
| `harmonized_rfu_qtl` | Unique normalized RFU-QTL records with evidence flags, target lists, counts, and tiers |
| `matched_evidence` | One RFU-QTL record × one external record per match; external fields have `evidence_` prefixes |
| `unmatched_variants` | Unmatched records from every layer, unresolved identities, and collapsed duplicate counts |
| `rfu_summary` | One row per `rfu_label`, evidence unions, gene/peak/trait lists, counts, maximum observed association tier |
| `variant_summary` | One row per build and ordered variant key, summarizing its RFU-QTL records |
| `credible_set_overlaps` | Pairwise descriptive RFU credible-set intersections with each external layer |
| `provenance` | Package/schema version, settings, row counts, resources, input metadata, scientific limitations |

`record_id` and `rfu_record_id` link pairwise evidence and QC back to normalized
records. IDs are deterministic SHA256 digests of normalized record content;
changing retained metadata changes identity. Identical normalized records collapse
with a `duplicate_record` QC count. Distinct effects, sources, targets, contexts,
or credible-set memberships remain distinct records. Conflicting estimates are
retained rather than pooled. There is no cross-layer Cartesian product.
Gene/peak/trait and shared-variant lists are sorted JSON arrays inside table cells,
so TSV exports retain unambiguous list boundaries. Table order is deterministic,
independent of input row order (not necessarily genomic position order).

Presence flags (`has_rfu_qtl`, `has_eqtl`, `has_caqtl`, `has_gwas`) describe supplied
records. `eqtl_finemapped`, `caqtl_finemapped`, and pairwise `evidence_finemapped`
mean a supplied PIP meets `pip_threshold` (default 0.95). Missing PIPs do not count
as support. Credible-set membership alone does not set these flags.
`shared_high_pip` requires the threshold on both RFU-QTL and external records.
These booleans indicate observed support, not whether a layer was adequately
measured; provenance records which inputs were supplied.

| Association tier | Evidence present |
|---|---|
| 1 | RFU-QTL, without eQTL/caQTL (possibly GWAS) |
| 2 | RFU-QTL and exactly one regulatory layer (possibly GWAS) |
| 3 | RFU-QTL + eQTL + caQTL |
| 4 | RFU-QTL + eQTL + caQTL + GWAS |

Tier 4 is **not causal**. `evidence_count` counts supplied layer types, including
RFU-QTL. GWAS-only support increases the count but leaves tier 1. Summaries report
`max_evidence_tier` over individual RFU-QTL records. An RFU with eQTL support at
one variant and caQTL support at another does not acquire tier 3. RFU/variant
summary flags are unions, so use pairwise records to inspect contexts and sources
before proposing a biological chain. Even within one association, layers may
originate in different conditions or populations.

No p-value filtering, multiple-testing correction, or significance discovery is
implicit. Supply evidence selected by an explicit upstream analysis, and record
its criteria in metadata. Feeding all nominal associations will count them all;
tiers are not strength-of-association scores.

## Credible-set overlap is separate from colocalization

```python
overlap = scrfu.tl.credible_set_overlap(
    rfu_qtl_members, eqtl_members,
    left_layer="rfu_qtl", right_layer="eqtl", pip_threshold=0.95,
)
```

A set is scoped by target (RFU/gene/peak/trait), source, release, context,
`locus_id`, build, and `credible_set_id`. IDs such as `cs1` must not ambiguously
identify multiple loci within that scope. Supply `locus_id` where IDs are reused.
Only set pairs with at least one **exact ordered variant-key intersection** are
returned. Unlike triangulation's optional allele reversal matching, reversed keys
are not intersected. No overlap means an empty output, not a negative statistical
colocalization result.

Outputs include supplied unique set sizes, shared members, Jaccard fraction, and
shared high-PIP members. Repeated membership rows count a variant once; if repeated
members have different PIPs, the maximum supplied PIP is used. Preserve distinct
fine-mapping runs using source/release/context/locus metadata. Counts refer only to
members supplied: do not prefilter complete credible sets to high-PIP variants if
you need meaningful full-set sizes. Record coverage/purity and upstream filtering;
set membership or a PIP does not establish a shared causal signal.

## Matos et al. local-file example

### Preparing the published RfuWAS RFU-QTL input

`examples/rfuwas_regulatory_prepare.py` converts the verified **Supplementary
Data 1** schema (`SNP`, `RFU`, `beta`, `t.stat`, `p.value`) into canonical
GRCh38 RFU-QTL evidence. It parses `CHR_POS_REF_ALT`, retains original columns,
preserves published one-based RFU identities as integer strings without an offset,
and records unresolved effect-allele orientation. It neither flips beta nor
reconstructs SE. Data 1 is not interchangeable with the **GRCh37 lasso prediction
weights**; **Data 4** contains RFU–phenotype associations to join by RFU rather
than variant. See the [verified resource audit](regulatory_reconnaissance.md).

```bash
python examples/rfuwas_regulatory_prepare.py \
  --input /path/to/42003_2024_7010_MOESM4_ESM.xlsx --format xlsx \
  --release s42003-024-07010-x --outdir /path/to/new/rfuwas_prepared
```

XLSX input explicitly selects `Data 1` and the headers on worksheet row 2. It
requires `openpyxl` in the analysis environment, not as a scRFU core dependency.
Alternatively, export **Data 1 only**, with its column header first, and use
`--format tsv` or `--format csv`; these need only pandas. Add `--nrows 25` for a
bounded smoke test. No files are downloaded or existing outputs overwritten.

Outputs are `rfu_qtl.tsv` and `provenance.json`, recording the source/release,
GRCh38 build, input SHA256, selected worksheet/export, row limit/counts, and
unresolved orientation. The complete prepared TSV can be used directly:

```python
rfu_qtl = pd.read_csv('/path/to/rfuwas_prepared/rfu_qtl.tsv', sep='\t',
                      dtype={'rfu_label': str})
result = scrfu.tl.regulatory_triangulation(
    rfu_qtl, eqtl=eqtl, caqtl=caqtl, allow_allele_reversal=False,
)
```

Remove `--nrows` when preparing the full Data 1 table for the post-reset run.
Direction comparisons remain missing until independently verified allele coding
is supplied. The callable example helpers are `prepare_rfuwas_data1` (DataFrame
conversion), `read_data1` (local reading), and `run` (reading plus export).

### Adapting Matos molecular-QTL files

The motivating [Matos et al. preprint](https://doi.org/10.64898/2026.01.27.26344979)
and its [public analysis repository](https://github.com/marlmatos/cd4t-qtl-map)
provide CD4 molecular-QTL mapping and fine-mapping workflows. A versioned
[data release](https://zenodo.org/records/18497815) is available separately. Obtain
permitted files yourself; scRFU does not download the archive. This example was
checked against the public scripts and synthetic exports, not the large released
archive. Actual release headers must be inspected before use.

`examples/matos_regulatory_triangulation.py` accepts a JSON manifest with local
TSV/CSV files (optionally gzip), source/release/build metadata, and explicit formats:

| Format | Expected columns / adaptation |
|---|---|
| `canonical` | Generic columns above; use for RFU-QTL/RfuWAS and GWAS |
| `tensorqtl` | `variant_id`, `phenotype_id` → gene or peak; optional `slope` → beta, `slope_se` → se, `pval_nominal` → pvalue, `af` → allele_frequency |
| `susie` | `gene` or `peak`; original `variant_id.ss` preferred over LD-aligned `variant_id`; `pip` or `PIP`, `cs` → credible_set_id, `region` → locus_id; `coverage` retained |

The study [ID-generation script](https://github.com/marlmatos/cd4t-qtl-map/blob/main/qtl_mapping/caqtl_mapping/005_Change_var_names.sh)
uses `CHR:POS[b38]REF,ALT`, while allele-parsing comments in the
[eQTL](https://github.com/marlmatos/cd4t-qtl-map/blob/main/qtl_finemapping_coloc/03.susie_finemap_eQTLs.R)
and [caQTL](https://github.com/marlmatos/cd4t-qtl-map/blob/main/qtl_finemapping_coloc/03.susie_finemap_caQTLs.R)
scripts describe a different order. Consequently custom IDs require explicit
`allele_order` (`ref-alt` or `alt-ref`), verified against the release VCF. No effect
allele is guessed from either order. Use a verified `effect_allele_column` or a
canonical `effect_allele` to enable direction comparisons. Embedded b37/b38 tags
must agree with the supplied build. Lifted coordinates must be provided as a
separate canonical input with its actual assembly, never by relabeling a build.

The fine-mapping scripts include a fallback with lower coverage. Preserve the
supplied `coverage`; do not assume all sets are 95% credible sets. `cs` is scoped
with `region`/target/context/source/release to avoid combining unrelated sets.

Example manifest (all paths relative to the manifest):

```json
{
  "rfu_qtl": {
    "source": "my-rfuQTL", "release": "analysis-v1", "genome_build": "GRCh38",
    "files": [{"path": "rfu_qtl.tsv", "format": "canonical"}]
  },
  "eqtl": {
    "source": "Matos-CD4", "release": "zenodo-18497815",
    "genome_build": "GRCh38", "context": "CD4-pseudobulk",
    "url": "https://zenodo.org/records/18497815",
    "files": [
      {"path": "eqtl.tsv", "format": "tensorqtl", "allele_order": "ref-alt"},
      {"path": "eqtl_credible_sets.tsv", "format": "susie", "allele_order": "ref-alt"}
    ]
  },
  "caqtl": {
    "source": "Matos-CD4", "release": "zenodo-18497815",
    "genome_build": "GRCh38", "context": "CD4-ATAC",
    "files": [{"path": "caqtl_credible_sets.tsv", "format": "susie", "allele_order": "ref-alt"}]
  }
}
```

`ref-alt` here illustrates a verified choice, not a default for every released
file. `files` may contain association and fine-mapping exports: they remain
separate evidence records so PIP-only rows are not assigned effects from an
unverified join. Row-level metadata is preserved; manifest source/release values
provide defaults. Per-file `sep: ","` supports CSV; Parquet/RDS exports must first
be converted locally to delimited tables (no new readers/dependencies are added).

For differing headers use a per-file `column_map` from source to canonical names,
plus `variant_column` when IDs have merge suffixes. Mapped columns are the names
used by `variant_column` and `effect_allele_column`. For example:
`"column_map": {"variant_id.y": "original_variant"}, "variant_column": "original_variant"`.
The adapter rejects conflicting coordinate IDs and unsupported ID forms.

```bash
python examples/matos_regulatory_triangulation.py \
  --manifest /path/to/local/manifest.json --outdir /path/to/new/results
```

Outputs are `regulatory_hits.tsv`, `rfu_regulatory_summary.tsv`,
`variant_regulatory_summary.tsv`, `unmatched_variants.tsv`, `provenance.json`,
plus `harmonized_rfu_qtl.tsv` and `credible_set_overlap.tsv`. Existing output files
are not overwritten. Each input gets a streaming SHA256 checksum; optional
manifest `sha256` values are verified. Paths, options, releases, builds, URLs,
checksums, input/unique row counts, and package/schema version are recorded.
There is no load timestamp, making repeated exports byte-reproducible for fixed
files, manifest, paths, software, and settings. Core DataFrame analyses retain
supplied resource/checksum metadata without claiming to verify unknown files.

## Integration and visualization

`join_regulatory_summary(existing_summary, result)` left-joins on `rfu_label`,
permits repeated RFUs on the left, preserves its order/index, and rejects column
collisions. Labels are matched as trimmed strings; original left labels remain
in output. Unobserved RFUs retain missing evidence. RFU summaries can therefore
join pseudobulk, phenotype coupling, longitudinal, antigen-evidence, and
frozen-reference transfer tables. Ensure RFU labels use the same frozen reference
and label namespace; string equality alone does not verify reference equivalence.

RFU-QTL/RfuWAS supplies the upstream variant–RFU association table, including its
model, covariates, sample size, ancestry, allele coding, abundance scale, and
multiple-testing decisions. This layer does not re-run RFU on CD4 expression data
or estimate the genetic association. Single-cell state, phenotype, and VDJdb
annotations provide additional orthogonal context. Regulatory overlap does not
prove antigen specificity, antigen recognition, or a causal phenotype mechanism.

Optional Matplotlib plots use the existing plotting extra:

```python
scrfu.pl.regulatory_evidence_heatmap(result)  # RFUs × observed evidence types
scrfu.pl.regulatory_evidence_bar(result)      # RFU counts per evidence type
```

The heatmap includes RFU-QTL, caQTL, eQTL, GWAS, and high regulatory PIP support.
These are RFU-level unions across variants/contexts, not a diagram of one shared
mechanism. Both functions accept `ax` and return an `Axes`; no network dependency
or automatic figure display is introduced.

## Future formal colocalization

A separate optional method would need region-wide, appropriately harmonized
summary statistics rather than a list of significant/shared variants: effect
estimates and variances (or method-supported alternatives), allele frequencies,
sample sizes, quantitative-trait variance or case/control proportions as needed,
and consistently defined loci. Multiple-signal models additionally require
ancestry-matched, allele-aligned LD, with validated matrix ordering and quality.
The method must specify priors, assumptions about causal-signal counts, population
and sample-overlap handling, and sensitivity analyses. Use a validated coloc or
SuSiE-coloc implementation behind an optional dependency and validate against
known examples and null simulations. Even formal colocalization is evidence for
a shared signal under a model; mediation and causal direction require additional
study design and assumptions. No formal colocalization method is implemented here.
