# Dataset validation plan

No accession is listed until its files, fields, license, and cohort design have
been verified. Dataset selection must not be driven by favorable results.

| Slot and role | Required fields and scale | Intended use / figure | Access and preprocessing | Confounding checks |
|---|---|---|---|---|
| Public bulk TCR aging cohort (development or validation) | donor, biological sample, age, chain, CDR3aa, V/J; enough donors for sample-level analysis | frozen-reference coverage, comparators, age-associated replication only after prespecification; Fig. 3 | accession/license unverified; AIRR/canonical conversion | age/sex/tissue/platform/depth/disease, repeated donors |
| Second independent or longitudinal TCR cohort | donor, sample, time, optional compartment, CDR3aa/V/J | transfer, longitudinal similarity/retrieval; Fig. 2 or 3 | unselected; explicit harmonization map | visit spacing, attrition, treatment, depth, batch |
| Paired single-cell TCR plus transcriptome/protein atlas | cell barcode, receptor chains, sample/donor, phenotype, optional modalities | phenotype linkage, scalability, bounded validation; Fig. 1/3 | Wells may fill this role; user-supplied H5AD/cache, no bundling | cell-quality, donor imbalance, first-N subset bias |
| Deep longitudinal validation cohort | repeated donor/sample/time and generic compartment | methods demonstration; Fig. 2 | separately governed runtime input; never bundled or used in tests | only six participants; no population-level aging claim |
| Completely held-out public TCR cohort | fields sufficient for frozen metrics, never used for tuning | final transfer/retrieval validation; Fig. 3 | must be registered and hashed before evaluation | cohort shift and missing harmonized fields |
| Optional public BCR cohort | paired heavy/light, V/J, isotype, SHM/germline, family, donor/sample | gated Figure 4 only | not selected; BCR gate must first pass | library chemistry, class switching, tissue, mutation calling |

The held-out role cannot be reassigned after results are viewed. Public licenses
must permit the intended redistribution or else only download instructions and
hashes may enter the reproducibility bundle.
