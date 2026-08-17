# BCR architecture design gate

## Current status

No BCR functional-unit model is implemented. The canonical receptor schema can
retain IGH, IGK, and IGL records, but applying a TCR-derived RFU reference to
those records is prohibited. This distinction is enforced in scientific scope,
even where generic adapters accept multiple chains.

## Minimum architecture

A receptor-specific BCR method must define and validate:

- separate IGH and IGK/IGL feature handling and productive-chain selection;
- explicit heavy/light pairing, including unpaired and multiply paired cells;
- V/J usage and CDR3 length representations appropriate to each chain;
- isotype and subclass, with missing or discordant constant-region evidence;
- somatic-hypermutation measurements with germline/reference provenance;
- naïve versus antigen-experienced and switched versus unswitched states;
- clonal-family identities and lineage structure without equating exact CDR3
  identity with a family;
- a receptor-specific reference construction and assignment algorithm;
- uncertainty and reference-coverage diagnostics;
- original/independent technical fidelity and held-out biological validation.

## Go/no-go gate

BCR is omitted from the first methods paper unless all conditions below are met:

1. a receptor-specific BCR reference and eligibility rule are frozen;
2. technical fidelity and processing-order invariance are demonstrated;
3. at least one suitable public BCR cohort is analyzed;
4. isotype or maturation-state associations are interpretable and not artifacts
   of chain selection or sequencing depth;
5. a completely held-out validation exists.

Until then, BCR is a data-model and architecture design topic only. A generic
sequence clustering placeholder would not satisfy this gate.
