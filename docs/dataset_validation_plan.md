# Dataset validation plan

No accession is listed until its files, fields, license, and cohort design have
been verified. Dataset selection must not be driven by favorable results.

| Slot and role | Required fields and scale | Intended use / figure | Access and preprocessing | Confounding checks |
|---|---|---|---|---|
| Public bulk TCR aging cohort (development) | donor, biological sample, age, chain, CDR3aa, V/J; enough donors for sample-level analysis | frozen-reference coverage, comparators, age-associated replication only after prespecification; Fig. 3 | PRJNA602091 selected; raw SRA scale and UMI-aware preprocessing manifest remain unresolved | age/sex/tissue/platform/depth/disease, repeated donors |
| Second independent or longitudinal TCR cohort | donor, sample, time, optional compartment, CDR3aa/V/J | transfer, longitudinal similarity/retrieval; Fig. 2 or 3 | GSE190905 receptor-only processed tables evaluated with explicit Pre/Post harmonization | treatment, depth, batch; six-patient technical cohort is not aging evidence |
| Paired single-cell TCR plus transcriptome/protein atlas | cell barcode, receptor chains, sample/donor, phenotype, optional modalities | phenotype linkage, scalability, bounded validation; Fig. 1/3 | Wells 1k/10k/25k deterministic subsets evaluated without loading expression matrices | cell-quality, donor imbalance, bounded samples are not claimed representative |
| Deep longitudinal validation cohort | repeated donor/sample/time and generic compartment | methods demonstration; Fig. 2 | separately governed runtime input; never bundled or used in tests | only six participants; no population-level aging claim |
| Completely held-out public TCR cohort | fields sufficient for frozen metrics, never used for tuning | final transfer validation; Fig. 3 | GSE157007 registered and all 17 VDJ input hashes frozen before RFU evaluation | cross-sectional one-sample-per-donor design prevents donor retrieval; three GEO assay labels conflict with their VDJ titles/files |
| Optional public BCR cohort | paired heavy/light, V/J, isotype, SHM/germline, family, donor/sample | gated Figure 4 only | not selected; BCR gate must first pass | library chemistry, class switching, tissue, mutation calling |

The held-out role cannot be reassigned after results are viewed. Public licenses
must permit the intended redistribution or else only download instructions and
hashes may enter the reproducibility bundle.

Verified accessions, publications, file availability, and role decisions are in
[`public_dataset_candidates.md`](public_dataset_candidates.md). Acquisition
details are frozen in [`acquisition_manifests/`](acquisition_manifests/).
