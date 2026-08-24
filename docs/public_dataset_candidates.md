# Verified public dataset candidates

Verified on 2026-08-24. Selection was based on receptor files and cohort design,
not observed scRFU results. NCBI places no restrictions on GEO molecular-data
reuse, but cannot transfer or adjudicate submitter intellectual-property rights;
scRFU therefore records download instructions and hashes rather than
redistributing source data.

## Ranked roles

| Rank / role | Accession and publication | Design | Receptor evidence | Access and restrictions | Decision |
|---|---|---|---|---|---|
| 1 — development aging cohort | PRJNA602091; *Longitudinal analysis reveals age-related changes in the T cell receptor repertoire of human T cell subsets*, JCI 2022, DOI 10.1172/JCI158122 | 30 healthy adults, two visits averaging 9.2 years apart; late 20s to early 80s at first visit; equal numbers of women and men; sorted naive/memory CD4/CD8 | TCRA and TCRB UMI-based RNA sequencing; raw sequence data are in SRA; publication supplies study metadata and analysis scripts | Public SRA. The raw-data byte size and run manifest were not frozen, and processed clonotype tables were not verified, so acquisition is blocked pending size and processing review. | Selected as development aging cohort; not acquired or evaluated. |
| 2 — independent technical validation | GSE190905 / PRJNA788905; PMID 39891774 | Six NSCLC patients, paired pre/post SABR or anti-PD-1 plus SABR; 43,051 metadata cells and 27,655 TCR-bearing cells | Processed wide TCR table contains primary/secondary TRA/TRB amino-acid CDR3, V/D/J, nucleotide junction, clonotype, and abundance fields | Public GEO; processed TCR table is 1.6 MB and metadata is 204.6 KB. No source data are redistributed. | Acquired and evaluated for frozen-reference transfer and donor retrieval; not an aging or population-level cohort. |
| 3 — held-out aging/frailty cohort | GSE157007 / PRJNA659762 / SRP279088; *Nature Aging* 2022, DOI 10.1038/s43587-022-00198-9 | 17 donors: 3 cord blood, 3 healthy young, 6 healthy old, 5 frail; 114,467 immune cells; cross-sectional | Single-cell TCR V(D)J. GEO records identify 17 VDJ samples and per-sample `filtered_contig_annotations.csv.gz` files produced by Cell Ranger 3.1.0 | Public GEO/SRA; full processed archive is 706.4 MB. The 17 receptor-only files total 4.5 MB and are recorded by hash without redistribution. | Preregistered as completely held out before receptor-outcome inspection; the immutable input manifest is frozen. |
| 4 — reserve validation cohort | GSE158848 / PRJNA666692 / SRP285949; PMID 33289628 | 86 repertoires from eight effector/memory CD4 subsets in five healthy donors; TRA and TRB, with selected biological replicates | Processed per-sample VDJtools-format clonotype tables; representative TRB file is 2.0 MB and includes MiXCR-derived assignments | Public GEO/SRA; combined processed TAR is 61.1 MB | Suitable reserve for phenotype/subset transfer, but too few donors for population claims and not an aging cohort. |

## Verification notes

- GSE190905 provides 16 sequencing records and four series-level processed
  files. Only the receptor and metadata tables were downloaded. The 55.1 MB UMI
  and 425.0 MB normalized expression tables were not downloaded.
- GSE157007 provides 48 sequencing records. The series metadata resolve to 17
  VDJ libraries. Fourteen are labelled as single-cell TCR; three later old-donor
  VDJ records are labelled as single-cell RNA in the assay field but have VDJ
  titles and `filtered_contig_annotations.csv.gz` files. This inconsistency is
  retained as a source caveat rather than silently corrected.
- GSE158848 explicitly separates TRA and TRB libraries. scRFU would use only
  TRB with the same official frozen reference.
- PRJNA602091 is scientifically strongest for aging/longitudinal validation,
  but raw-data scale and the absence of a verified processed clonotype bundle
  prevent automatic acquisition in this session.

## Authoritative sources

- [PRJNA602091 publication and data availability](https://pmc.ncbi.nlm.nih.gov/articles/PMC9433102/)
- [GSE190905 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190905)
- [GSE157007 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157007)
- [Representative GSE157007 VDJ sample](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4750312)
- [GSE158848 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158848)
- [Representative GSE158848 TRB sample](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4812080)
- [NCBI GEO data disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html)
