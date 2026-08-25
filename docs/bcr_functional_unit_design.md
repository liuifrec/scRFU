# Experimental BCR functional-unit design

Status: design only; no validated or production BCR functional-unit reference
exists. Canonical TCR RFU assignment remains TRB-only and its frozen definition
is unchanged.

## What a frozen BCR representation must preserve

A transferable BCR state should retain receptor recognition constraints while
not merely rediscovering donor-specific clonal lineages. Its definition must be
deterministic, independent of phenotype/outcome labels, explicit about missing
light chains and mutation fields, and stable under donor and cohort holdout.

| Representation | Biological interpretation | Required data | Transfer and missing-data behavior | Main bias / limitation |
|---|---|---|---|---|
| A. Sequence-only grouping | Similar junction sequences, predominantly recognition-related but also lineage-related | CDR3 amino-acid sequence and chain | Broadest availability; heavy and light chains must be separated; missing light is irrelevant in heavy-only mode | Can merge unrelated receptors or split matured relatives; ignores gene context and pairing |
| B. V/J plus CDR3 similarity | Gene-contextualized junction similarity | Chain, V, J, CDR3 | Portable when gene calls are normalized; missing V/J requires an explicit lower-information stratum, never imputation | Annotation-tool and allele-resolution bias; still conflates recognition and lineage |
| C. Heavy-chain receptor state | IGH recognition and rearrangement state | Productive IGH with V/J/CDR3 | Works for bulk and single-cell cohorts; light-chain absence is explicit | Cannot represent light-chain contributions and may overstate similarity |
| D. Paired heavy/light grouping | More complete receptor recognition state | Natively paired productive IGH plus IGK/IGL | Highest specificity for paired single-cell data; heavy-only records require a distinct coverage class, not forced pairing | Poorer cross-study coverage and sensitivity to multiplets/chain selection |
| E. Maturation-aware grouping | Receptor state plus class switching and affinity maturation | Above receptor fields plus explicit isotype/constant region and SHM | Missing isotype/SHM must remain unknown; sensitivity analyses must compare receptor-only and maturation-aware variants | Isotype, tissue, assay and SHM can dominate distance and reduce portability |
| F. Clonal-family-aware representation | Lineage maturation and diversification | Nucleotide sequence, germline annotation and a prespecified family definition | Valuable within donor/time series; family identifiers are cohort-local and cannot be assumed transferable | Clonal families are not functional units and are strongly donor-specific |

## Proposed evaluation sequence

1. Validate the canonical BCR table and field provenance on a compact public
   10x/AIRR dataset.
2. Measure heavy/light pairing, isotype, SHM and family-field coverage without
   using biological outcomes.
3. Prespecify and compare B and C as lower-information baselines and D as the
   paired model. Treat E as a sensitivity layer, not the default distance.
4. Keep F as a lineage comparator and longitudinal summary; never rename a
   clonal family as a functional unit.
5. Build any prototype on training donors only and freeze its source hashes,
   normalization, distance, seed and coverage rules before donor/cohort holdout.
6. Require deterministic rebuild, missing-feature sensitivity, downsampling
   stability and independent public-cohort transfer before exposing a stable
   API.

## Missing-data policy

- Heavy-only, paired-heavy/light and light-only records are distinct coverage
  states.
- Missing V/J, isotype or SHM is represented as unavailable, not a biological
  category and not silently imputed.
- Class switching derives only from explicit constant-region/isotype evidence.
- Multiple productive chains are retained with deterministic selection
  provenance; they are not silently discarded.
- Reference coverage must be reported separately for receptor-only and
  maturation-aware models.

## Go/no-go gate on 2026-08-25

**NO-GO for BCR reference construction.** Public candidates have been verified,
and experimental preprocessing/state-feature code exists, but no candidate has
yet completed local field-completeness QC and no outcome-free distance/reference
construction has been prespecified and independently held out. Implementing a
reference now would violate the required gate and risk turning clonal or assay
structure into an alleged transferable unit.

The package may continue to ship clearly experimental BCR preprocessing and
interpretable state summaries. A later prototype can proceed only after the
first three evaluation steps above are satisfied.
