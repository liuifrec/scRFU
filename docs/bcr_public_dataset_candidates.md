# Experimental BCR public-data candidates

Verified on 2026-08-25 from the linked public repository records. This is a
software-development shortlist, not a claim that scRFU has a validated BCR
functional-unit reference. TCR RFU centroids must never be applied to these
data.

## Ranked candidates

| Rank / use | Public record | Design and size | Verified receptor content | Isotype / SHM / family information | Access and caveats |
|---|---|---|---|---|---|
| 1 — longitudinal method development | GSE219098 / PRJNA907025 / SRP410707; PMID 37440409 | Three vaccinees sampled longitudinally through six months; spike-binding and non-binding B cells plus plasmablasts; 81 GEX/ADT/VDJ records | Single-cell 10x BCR libraries with paired receptor reconstruction; processed CSV/FASTA archive is 157.4 MB | The study reports progressive SHM and tracked clones; phenotype comes from sorting and paired single-cell measurements. Exact field completeness must be measured after acquisition. | Public GEO/SRA. Strongest development design, but donor count is small and antigen enrichment means it is not a general-population reference. |
| 2 — maturation/phenotype validation | GSE244297 / PRJNA1022135; PMID 39168129 | Five volunteers sampled before and at weeks 1–3 after two vaccinations; one donor excluded by the original study; 24 CSP/GEX/VDJ records | Eight single-cell 10x VDJ libraries in a 117.7 MB processed archive | Study design distinguishes spike-reactive, activated B cells and plasmablasts and analyzes affinity maturation. Exact isotype, SHM and family-column coverage requires local QC. | Public GEO/SRA. Sorting and antigen enrichment are useful validation labels but must not enter reference construction. |
| 3 — compact adapter and pairing validation | GSE266519 / PRJNA1107496; PMID 40294012 | Four transplant recipients and three healthy controls with pre/post-vaccination and some post-infection samples; 48 modality records | One 3.0 MB `filtered_contig_annotations.csv.gz` BCR table and one 535.7 KB clonotype table from paired 10x VDJ libraries | 10x annotations supply heavy/light contigs, V/D/J/C calls, productivity and clonotype IDs; SHM and clonal-family fields are not promised by the repository record | Public GEO/SRA. Compact processed files make this the preferred real adapter fixture, but the disease design is unsuitable for an outcome-free general reference. |
| 4 — redistributable tutorial/adapter check | 10x Genomics human B cells pre/post influenza vaccination v2 | One healthy donor, pre-vaccination and day 14; 18,115 cells | Paired BCR and gene-expression libraries analyzed with Cell Ranger 5.0.0 | Constant-region and clonotype fields are expected in Cell Ranger VDJ output; SHM and inferred clonal families are not advertised | CC BY 4.0. Useful for format testing only; one donor cannot support transfer validation or reference construction. |
| 5 — broad sequence resource | Observed Antibody Space (OAS), paired collection | Paired sequences from multiple public 10x studies with subject, disease, vaccination, tissue, B-cell type and longitudinal metadata where supplied | Paired heavy/light AIRR-like annotated CSV files with V/D/J calls, junctions and germline alignments | Germline alignments permit mutation analyses; metadata completeness varies by contributing study | OAS states CC BY 4.0. Cohort boundaries and source-study hashes must be frozen before use; the whole collection must not be downloaded merely for exploratory prototyping. |

## Field and suitability conclusions

- GSE219098 is the highest-priority longitudinal BCR development cohort because
  the public record verifies repeated time points, donor demultiplexing, paired
  single-cell BCR and an SHM-focused design.
- GSE244297 is the strongest phenotype/maturation validation candidate, but its
  sorted biological labels cannot be used to build an outcome-free reference.
- GSE266519 is the safest first real-data adapter test because the processed BCR
  contig table is small. It does not by itself justify a frozen BCR reference.
- OAS is a possible future multi-study source only after explicit cohort
  selection. Its scale and heterogeneous metadata make indiscriminate download
  inappropriate.
- No candidate is approved here as a BCR functional-unit reference. Approval
  requires field-completeness QC, a prespecified receptor distance, donor-aware
  construction, defined missing-data behavior and an independent holdout.

## Completed bounded real-data QC

The processed GSE219098 and compact GSE266519 receptor tables were acquired
outside Git. The GSE266519 gzip SHA256 is
`ce615c2ca41ed460617328f1b0280d8ad348974e6636eb156e6174cd7a562f3c`.
It contains 34,477 productive contigs from 16,143 cells: 16,852 IGH, 10,636
IGK and 6,989 IGL rows. Deterministic pairing yielded 15,494 paired, 177
heavy-only and 472 light-only cells. Explicit/constant-derived isotype was
available for 16,837 rows; mutation-frequency and clonal-family fields were
absent. The file remains external and no source rows are redistributed. These
GSE219098 yielded 379,741 productive contigs from 177,663 cells, including
171,642 paired cells. Neither dataset supplied SHM/germline-identity or inferred
clonal-family fields in the acquired tables. Exact outcome-independent
representations transferred to only 0–0.236% of eligible GSE266519 cells. These
results validate preprocessing and feature extraction but reinforce the
BCR-reference NO-GO documented in `bcr_feasibility_report.md`.

## Authoritative sources

- [GSE219098 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE219098)
- [GSE244297 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244297)
- [GSE266519 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266519)
- [10x pre/post influenza B-cell dataset](https://www.10xgenomics.com/cn/datasets/human-b-cells-from-a-healthy-donor-pre-and-post-flu-vaccination-v-2-2-standard-5-0-0)
- [OAS paired-data documentation](https://opig.stats.ox.ac.uk/webapps/oas/documentation_paired/)
- [OAS licensing statement](https://opig.stats.ox.ac.uk/webapps/oas/)
- [AIRR Data Commons access documentation](https://docs.airr-community.org/en/latest/adc/data_submission.html)
