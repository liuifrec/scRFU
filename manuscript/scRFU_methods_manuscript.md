# scRFU: reference-anchored analysis of longitudinal T-cell repertoire remodeling with applications to radiotherapy

Manuscript v1, 2026-09-16. This is a results-bearing draft, not a submission-ready
manuscript. RP1-14 reconstruction and its six-panel longitudinal figure are
complete; the external GSE280982 endpoint application remains unfinished. The paused unpublished survivor cohort contributes
no observations, outcome-selected features, transformations or thresholds.

## Abstract

Longitudinal immune-repertoire comparisons depend on which receptors are
observed, how expanded clones are weighted and how sequences are grouped.
We describe scRFU, a reproducible interface to the published frozen-reference
Repertoire Functional Unit (RFU) method, with explicit assignment provenance
and longitudinal measurement at receptor, group and cellular-annotation scales.
After repairing a demonstrable row-alignment error in historical RP1-14 exports,
we reused assignments from six donors, three visits and both CD4/CD8 compartments
over 14.95–24.86 years. Median receptor/RFU total variation was 0.825/0.305 in
CD4 and 0.718/0.431 in CD8. Fixed matched groupings produced similar aggregation
effects. Using existing assignments from a public radiotherapy dataset, we also
recovered six paired donors and 27,655 TCR-bearing cells. On threshold-qualified cells, median
pre/post total variation was 0.840 for primary-TRB identities and 0.733 for RFUs;
median within-donor aggregation cancellation was 0.099. Group-size and receptor-
feature-matched controls produced similar cancellation, while depth and clone
weighting changed the measured distances. These findings demonstrate the need
to report aggregation and sampling effects, without attributing group stability
to preserved antigen-specific function. We also reused 303,088 primary-TRB
assignments from the Wells atlas for donor- and tissue-resolved cellular context,
and retained a completed regulatory-evidence example with explicit limits on
colocalization and causality. The contribution is a validated analysis workflow;
neither the original RFU method nor total-variation contraction is new here.

## Introduction

Exact receptor sequences provide a natural description of clonotypes, but sparse
sampling and clonal expansion complicate comparisons across visits and people.
Grouping receptors into a reference-defined space can support comparable
summaries while hiding changes among the group's constituent receptors. A
longitudinal method should therefore preserve the relationship between observed
clone change, group change, sampling coverage and independently annotated cell
state.

RFUs were introduced by Hu and colleagues in
[Quantifiable blood TCR repertoire components associate with immune aging](https://doi.org/10.1038/s41467-024-52522-z).
scRFU implements reproducible access to that external frozen reference within
pandas and scverse workflows. Its contribution is transparent data handling,
assignment fidelity, reusable longitudinal measurements and explicit evidence
boundaries. Reference membership alone does not prove a shared antigen target.

Radiotherapy provides a useful application because relative repertoire change
can reflect receptor redistribution, state composition, expansion and sampling.
Here we use published or public data only, with separate roles for a
published RP1-14 longitudinal demonstration, an existing multi-tissue atlas and
public radiotherapy studies. The present draft reports completed public-data
measurements and marks unavailable analyses instead of substituting findings
from an unpublished cohort.

## Results

### Existing assignments support a reproducible public-only analysis

The recovered GSE190905 and Wells runs both used the standard upstream
`AssignRFUs` backend, the same frozen reference and threshold 0.6. Their saved
receptor outputs passed completed-run validation and recorded checksum checks.
We did not rerun receptor assignment. RFU labels retain the original `RFU` prefix
and one-based numbering. The older Wells map-aware cache was unavailable and
was not treated as equivalent merely because its atlas filename was similar.

### Multiscale change distinguishes receptor turnover from aggregation

We apply one fixed receptor-to-group map at every visit, normalize each observed
sample consistently and compute total variation (TV) before and after grouping.
The difference, aggregation cancellation, is change hidden by the grouping.
Comparison with V/J and constrained random maps separates this mathematical
effect from claims about RFU-specific biology. Missing visits remain missing;
coverage and clone/read weighting define the estimand, rather than being hidden
inside a group-level score.

### RP1-14 recovers long-term donor profiles despite constituent-receptor turnover

The restored published-only material contains six donors, three chronological
visits and both CD4/CD8 compartments at every visit: 36 repertoires, 4,342,203
sequence rows and 70,960,689 sequencing reads. Collection-age metadata establish
14.95–24.86 years of follow-up. Source filename suffixes run in reverse time
order; they were not used as chronological labels. These are the healthy-volunteer
longitudinal samples described by
[Yoshida et al.](https://doi.org/10.1016/j.exger.2017.05.015), not the paused
144-person survivor cohort or a radiotherapy exposure comparison.

The historical 5,000-by-36 RFU matrix counts assignment rows per 10,000, including
below-threshold nearest labels. Its accompanying receptor exports contained a
row-alignment error: source labels were truncated while the upstream encoder
filtered CDR3s not starting with C. Exact source reconciliation repaired 330,433
row labels across 359,530 records without changing a single RFU or score. Repeated
sequence assignment conflicts disappeared, and all 216 receptors in a bounded
standard-reference check reproduced their historical labels. The recovered
map-aware function contains this truncation; its original invocation manifest
is unavailable. The [reuse audit](../docs/rp1_14_reuse_audit.md) preserves the
checks and remaining provenance limit.

A fixed dictionary of historical sequence assignments was applied to each
visit's full productive counts. Thus a clone falling below the old top-10,000
cutoff was not automatically treated as absent. The primary threshold-qualified
analysis covers 508,972 distinct nucleotide-CDR3/V/J receptor identities,
617,506 sample–receptor observations, 4,958 RFUs and 25,513,582 reads. Per-sample
coverage is 28.35–78.44% of productive reads, or 23.55–58.09% of all source reads.
Results remain conditional on the recoverable assignment universe. They do not
represent a new complete-repertoire assignment run.

| Compartment | Donors | Median receptor TV | Median RFU TV | Median cancellation | Median V / V–J TV |
| --- | ---: | ---: | ---: | ---: | ---: |
| CD4 | 6 | 0.8246 | 0.3055 | 0.4640 | 0.0549 / 0.1229 |
| CD8 | 6 | 0.7180 | 0.4309 | 0.2867 | 0.2034 / 0.2864 |

These earliest-to-latest comparisons use the same receptors for every grouping.
All 432 repeated donor/visit-pair/policy/weighting/grouping rows satisfy TV
contraction. CD8 RFU TV exceeds paired CD4 TV in five of six donors; the median
paired difference is +0.1192. This describes the observed small cohort, not a
population aging model. Different grouping granularities prevent interpreting
smaller V/J distances as greater functional stability.

Across all observed visit pairs, equal-donor mean within-person RFU cosine
similarity is 0.7609 in CD4 and 0.6763 in CD8, versus 0.3033 and 0.0426 for each
donor's average between-person comparisons. At a five-read detection threshold,
median persistent/union RFU fractions are 0.7942 and 0.7400. Among persistent
RFUs, median fractions with no shared observed nucleotide/V/J receptor are
0.4140 and 0.5983. Persistent groups can therefore contain different observed
receptors; neither detection nor persistence establishes antigen function.

Expansion and depth strongly affect the measurement. Unique-receptor weighting
reduces median RFU TV to 0.2365/0.2922 (CD4/CD8). Removing each pair's dominant
receptors gives 0.3053/0.4117. At the deliberately stringent common depth of
500 reads, median within-donor subsample RFU TV is 0.8235/0.7595. Those 50
read-subsampling replicates are neither cell/template resampling nor donor
confidence intervals. Thirty fixed feature-matched groupings give median
cancellation 0.4786/0.2942, close to the observed 0.4640/0.2867; this does not
support special functional stability of RFUs. Legacy age-related RFU scores and
mixed-cohort embeddings were not reused.

### GSE190905 measures radiotherapy-associated remodeling with sampling dependence

The cached GSE190905 release contains 43,051 RNA metadata cells and 27,655
TCR-bearing cells from six donors with two visits each. Four donors have the
source treatment label `SBRT` and two `I-SBRT`. Every TCR cell joined uniquely
to source RNA metadata and the saved receptor assignments. All donor, visit and
treatment labels agreed across the paired tables. Eight physical library pools
were reconciled to GEO records by exact annotated donor/visit membership;
processed barcode-prefix numbering differs from GEO Batch numbering. We did
not infer patient identity from shared receptor sequences or treat a pooled
library as one patient.

The [final publication](https://doi.org/10.1007/s00262-024-03935-8) describes
seven patients and 57,738 cells, blood collected before and after SABR, and a
prior immunotherapy-plus-chemotherapy group. Our results apply to the cached
six-donor release; the absent seventh donor is not reconstructed. GEO's
metastatic-disease wording also differs from the paper's early-stage description.
These discrepancies limit clinical subgroup interpretation.

The primary analysis retained 21,657 threshold-qualified cells. Assignment
coverage among TCR-bearing cells ranged from 72.45% to 91.82% across the twelve
samples. All six donors met the prespecified minimum of 100 retained cells per
visit for total, CD4 and CD8 comparisons. Source expression clusters had variable
support: for example, CD8 effector comparisons were available in five donors,
whereas the CD4-naive comparison had only two. Unsupported strata were
recorded without imputing visits.

| Grouping | Donors | Median receptor TV | Median group TV | Median aggregation cancellation |
| --- | ---: | ---: | ---: | ---: |
| RFU | 6 | 0.8399 | 0.7332 | 0.0986 |
| TRBV | 6 | 0.8399 | 0.2419 | 0.5675 |
| TRBV/TRBJ | 6 | 0.8399 | 0.4524 | 0.3611 |

Cancellation was computed within each donor before summarizing; the median
difference need not equal the difference of medians. The groupings are not
nested and have different resolutions. Smaller distances for coarse V/J
groupings do not demonstrate greater biological stability. The contraction
inequality held for all 462 reported donor/subset/policy/weighting/grouping rows.
Those rows represent repeated analyses of six donors, not 462 participants.

Threshold-qualified RFU cosine similarity averaged 0.4997 within donor and
0.0547 for each donor's mean between-donor comparisons. The pattern was
heterogeneous: two donors had very low within-person similarity. Shared
between-donor comparisons are dependent; no pair-level significance test was
used.

Across 30 fixed random maps matched on group size, V, J and CDR3-length bin,
98.82% of eligible unique receptors lay in strata able to exchange labels.
Approximately 90.0–91.1% of labels changed. Observed RFU cancellation was close
to the control distributions for several donors, with departures in both
directions. These descriptive controls do not establish a distinctive functional
stability property of RFUs or a calibrated biological null.

Empirical subsampling of actual cells to equal pre/post depth, capped at 500
cells, generally increased observed RFU distances. One donor had only 404
assigned cells at the smaller visit. Fifty replicate ranges quantify instability
under observed-cell subsampling, not uncertainty across participants. With
unique-receptor weighting, median total-repertoire RFU TV was 0.8215 compared
with 0.7332 under cell weighting. Removing the union of each visit's dominant
primary-TRB identity retained only 36.9% of one donor's post-treatment assigned
cells, illustrating substantial expansion sensitivity.

Within CD4 and CD8 categories, median cell-weighted RFU TV was 0.8166 and 0.6702,
respectively; median cancellation was 0.1187 and 0.0120. These are descriptive
paired measurements with different depths, expansion patterns and receptor
diversity, not evidence that radiation has a larger causal effect on CD4 cells.
Observed RFU persistence, constituent-receptor sharing and dominant-receptor
fractions are retained in the source tables at both one-cell and five-cell
detection thresholds. Nondetection remains sampling-dependent.

### Wells provides cellular context without selecting survivor-associated RFUs

The recovered atlas contains 610,429 cells from 24 donors. Its completed standard
assignment run contains 303,088 primary-TRB cells from 21 donors across ten
tissues; 233,913 cells pass threshold and represent 4,928 RFUs. The three atlas
donors without eligible primary TRB remain in the coverage denominator. We
saved donor/tissue/RFU/cell-type counts and within-group observed receptor
dominance for all eligible RFUs. A public-data support rule of at least 50 cells
across four donors retained 1,471 RFUs, without reference to survivor outcomes.

This supports linking a fixed receptor representation to source cellular
annotations while retaining donor and tissue denominators. The source
`cell_type` field supplies the main cellular categories; `cell_state` is mostly
unannotated and is not interpreted as a complete state taxonomy. Cell counts
are not combined with bulk read/template abundances. CDR3 amino-acid plus V
identity in these tables is an observed receptor proxy, not a verified paired-
chain or nucleotide-defined clonal lineage.

### Frozen supporting evidence does not establish prediction or mediation

The completed regulatory application is reused unchanged from `d62f1c1`.
Among 623 published RFU-QTL variants, 17 overlap Matos conditional eQTL
selections and 32 overlap conditional caQTL selections. Six **unique variants**
overlap both conditional layers at the same variant. They connect to 62
variant–RFU pairs and 36 RFUs, with six eQTL and ten caQTL target/rank records;
they are not six statistically independent mechanisms. Molecular association,
conditional selection, fine mapping and colocalization remain distinct, and
formal RFU colocalization is not established.

The previously completed GSE190905 prediction benchmark also remains unchanged.
Across six held-out donors, clone-weighted mean log loss increased from
1.850996 to 1.855212 when RFU categories were added to the TRBV/TRBJ/CDR3-length
baseline (delta +0.004216). A stricter analysis excluding receptor identities
shared across training/test donors gave 1.852679 versus 1.856996 (delta +0.004317).
No improvement was demonstrated in either fixed design. This does not invalidate
longitudinal measurement utility, nor does it establish universal absence of
RFU–state relationships. No classifiers were refit for this manuscript.

No radiotherapy–regulatory shortlist intersection is claimed here. Linking RFU
labels across a published genetic analysis and the reused single-cell reference
still requires reference-identity confirmation. Published RfuWAS disease links
remain downstream genetic-prediction annotations, not independent experimental
validation or evidence that molecular QTLs mediate radiation effects.

### Explicitly unfinished applications

**GSE280982:** public metadata list paired processed GEX/TCR resources for eight
tumor visits across three donors, plus three paired blood visits. The first
tumor donor has GEX entries without released matching TCR entries in the
inspected record, and one later tumor visit is absent. These are resource-
availability findings, not external endpoint results. The
[source publication](https://doi.org/10.1038/s41467-025-60827-w) concerns tumor
biopsies around radiotherapy; blood and tumor analyses will be reported as
portability applications, not interchangeable effect replications. No external
RFU assignments have been started in this checkpoint.

## Methods

### Frozen assignment and analysis boundary

We validated cached assignment tables and manifests rather than recomputing
RFUs. The common `km5000noMax.Rdata` SHA256 is
`64783074360edeca84b3f49ab9d682263e954752fc09bb208edaabde0fd82553`.
The backend, RFU code, trimer reference, threshold and source-output hashes are
recorded in the configuration and external provenance. Completed-run validation
checks chunk completion, reconstruction, labels and source checksums.

Only allowlisted datasets with explicit authorized sample membership enter the
manuscript runner. Source files are pinned by checksum. The unpublished
144-person cohort and mixed-cohort coordinates, scores and outcome-based
shortlists are excluded. Restored RP1-14 source files are protected by exact local
Git-ignore paths; they are not staged or redistributed. All source crosswalks and
derived biological tables remain outside Git.

### RP1-14 source reconciliation and count units

A study-specific loader verifies metadata/filename compartment agreement,
chronology by collection age, complete source sample membership and both read-
count header variants. Source percent frequencies must reconcile to the complete
integer-read denominator. CDR3 nucleotide intervals are sliced using the
zero-based source V index and CDR3 length, then translated and checked against
the supplied amino acids. Receptor identity uses that nucleotide CDR3 and source
maximum-resolved V/J calls. Family-only/unresolved gene categories are retained
explicitly in the V/J grouping controls.

Historical sequence/V prefixes, upstream C-start filter lengths, threshold
flags, matrix reconstruction and repeated-sequence consistency gate alignment
repair. The existing RFU vectors remain unchanged. A deterministic six-receptor-
per-sample comparison checks bounded current-reference parity; it is not a new
assignment campaign. The recovered union of productive historical amino-acid
assignments is applied unchanged to the full productive counts at all visits,
with coverage against both productive and all-source read denominators. Unmapped
mass is never combined into a fictitiously stable RFU.

Primary comparisons use threshold 0.6, earliest/latest visits and source read
weights. All observed visit pairs, nearest labels and unique-receptor weights
are sensitivities. The frozen development control settings are retained: 30
fixed maps permuted within V/J/five-amino-acid-length bins, 50 observed-count
subsamples capped at 500, and dominant-receptor removal. For this bulk dataset,
hypergeometric subsampling draws actual reads without replacement; it neither
estimates independent template sampling nor removes PCR dependence. The six
people, rather than reads, RFUs, visit pairs or control maps, are the biological
replicates. The RP-specific adaptation was recorded before endpoint inspection.

### Longitudinal design and grouping

GSE190905 cells were joined by unique source barcode. Donor and visit labels
were checked across the source TCR and RNA metadata. Pooled libraries were
reconciled against GEO sample titles by their full donor/visit membership, not
by similarly numbered prefixes. This reuses the authors' sample demultiplexing;
raw HTO data were not independently reprocessed.

The primary observed receptor identity is primary-TRB nucleotide CDR3, source V
and source J. All comparisons use a fixed mapping across both visits. Relative
mass is normalized separately within each observed sample. Primary analyses
condition on threshold-qualified cells and report excluded mass; nearest-policy
comparisons are separate sensitivities. Zero-mass samples have undefined
distances. Missing visits are never silently imputed.

We calculate total variation before and after grouping and report their
nonnegative difference as aggregation cancellation. This is the standard
contraction of total variation under a deterministic map. We make no claim of
mathematical novelty. Synthetic tests cover sampling from unchanged latent
frequencies, one-clone expansion, within-group replacement, between-group
redistribution, changing coverage, identical/empty observations and missing
visits. Biological replication is at the donor level.

The [frozen analysis specification](radiation_methods_specification.md) gives
the exact support thresholds, feature-matching constraints, seeds, depth and
dominant-receptor sensitivities. Subsampling uses actual cells without
replacement. Synthetic multinomial sampling is identified as simulation, and
integer counts are never inferred from transformed abundances. Source expression
categories define cell states independently of RFU labels; RFU labels are not
used to relabel cells. The earlier prediction benchmark's limitations concerning
receptor-gene contributions to expression labeling remain in its evidence record.

### Reproducibility and completion

The runner saves source tables, figures, software versions, executable/import
paths, input/configuration/code hashes and the command. A completion manifest is
written only after all scoped outputs exist. Reuse verifies both the scientific
fingerprint and every listed output hash. Preparation, measurement and export completion are recorded separately so an
interrupted figure/summary cannot masquerade as a finished analysis. Completed
RP1-14 measurements are reused when only exports require repair. Unexecuted
GSE280982 endpoints remain distinct from completed development results.
The original failed library-number assumption was corrected by source membership
reconciliation. A coverage-denominator test preserves Wells donors without
eligible TRB; corrected version `development_v1_1` leaves all GSE190905 source
tables byte-identical to the initial completed run.

## Discussion and limitations

The completed application demonstrates reproducible measurement at different
repertoire scales, not RFU superiority. Grouping necessarily hides some changes;
its extent depends on grouping resolution and the observed receptor universe.
The matched-group controls and depth sensitivities argue against interpreting
small group distances as preservation of antigen-specific immunity. Changing
assignment coverage and dominant-clone weighting can also alter longitudinal
summaries.

The radiotherapy results are observational relative-composition measurements
from six donors. Prior systemic therapy, sample depth, clonal expansion and
source-release differences preclude isolating a causal radiation effect. They
do not measure absolute lymphocyte depletion. The development dataset was
previously examined, and no external endpoint replication is claimed. Wells is
cross-sectional and contains heterogeneous tissue and donor support, not a
longitudinal treatment replication.

RP1-14 is a six-person longitudinal benchmark, not radiation biology or a
population aging model. Its recoverable assignment coverage is incomplete,
source reads can be PCR-dependent, and historical-reference identity is supported
by a bounded parity check rather than an original execution manifest. Formal RFU colocalization requires
dense regional RFU-QTL statistics, verified effect alleles and suitable study-
specific LD; significant-only associations and GRCh37 lasso weights cannot
substitute. No antigen, functional-recovery, genetic-mediation or clinical-
prediction claim follows from the present evidence.

## Data and code availability

Code and synthetic tests are on `manuscript/radiation-methods`. Public sources
are [GSE190905](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190905),
[GSE280982](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280982) and the
previously obtained [Wells atlas](https://doi.org/10.1038/s41590-025-02241-4).
The H5AD's own citation identifies the T-cell dataset version
`965008a7-e698-413c-a746-8855a945a7c5`, titled “Immune aging project - t cell
subset,” in [CELLxGENE collection cc431242](https://cellxgene.cziscience.com/collections/cc431242-35ea-41e1-a100-41e0dec2665b).
Its cached file hash and assignment provenance are recorded. The
external RFU reference and biological datasets retain their original access and
licensing conditions. No biological tables, donor crosswalks or reference assets
are redistributed in Git. RP1-14 access and redistribution are not presumed
from publication alone. The regulatory resources and limitations remain in the
[completed report](../docs/regulatory_realdata_pilot.md).

## Figure legends and source index

**Radiotherapy longitudinal figure** (`gse190905_longitudinal.pdf`). A, donor-
paired primary-TRB and group total variation on identical threshold-qualified
cells. Lines connect the same donor across groupings; the groupings are not a
nested biological hierarchy. B, threshold coverage among TCR cells before and
after treatment. C, one within-person cosine value and one mean between-person
value per donor. Points do not represent independent cell or sample-pair tests.

**Multiscale sensitivity figure** (`gse190905_multiscale_sensitivity.pdf`). A,
observed RFU distance (blue) and median/2.5th–97.5th percentiles of 50 empirical
equal-depth subsamples (gray). B, observed cancellation and the corresponding
range across 30 fixed size/feature-matched random groupings. These ranges are
not donor confidence intervals or calibrated null tests.

**RP1-14 longitudinal figure** (`rp1_14_longitudinal.pdf`). A, six donors with
three CD4/CD8 visits; marker area indicates threshold-reuse coverage among
productive reads. B, donor-level means of within- and between-person RFU cosine
similarity. C, receptor versus RFU TV on identical mapped productive reads.
D, paired CD4/CD8 earliest/latest RFU TV. E, persistent/union RFU detection and
the fraction of persistent RFUs without shared observed receptors (at least five
reads per visit). F, observed cancellation and 30 fixed feature-matched map
controls. Neither control ranges nor repeated sample pairs are independent
participants. Persistence is an observed grouping property, not antigen-function
preservation.

Separate external `figure_source_index.tsv` files map all five development panels
and six RP1-14 panels to source hashes, exact commands, denominators and limits.
Framework and external endpoint panels remain unassembled. The
[execution state](../docs/radiation_methods_execution_state.md) and
[claim-to-evidence table](radiation_methods_claims.tsv) identify the active
checkpoint and unfinished work.
