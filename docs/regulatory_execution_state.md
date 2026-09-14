# Regulatory manuscript application: execution state

Updated 2026-09-14. Branch `feature/regulatory-triangulation`; starting checkpoint
`ce2359e` is already on origin. Historical outputs and the pilot report are retained.
External workspace: `/home/liuyuchen/data/scRFU_regulatory_20260909/`.

## Analysis specification

1. **A — shared variants:** exact GRCh38 identities; no allele reversal, liftover,
   or inferred effect allele. Preserve marginal and conditional statistics.
   Query Bonferroni family is unique molecular variant–target tests before RFU /
   disease joins. Source marginal significance remains unknown without thresholds.
2. **B — shared regulatory signal:** source within-target conditional selections
   are identified by layer, context, target and rank. They are not independent
   mechanisms across targets. Source molecular coloc must identify its target /
   signal pair, posterior, membership and provenance; proximity is separate.
3. **C — RFU beyond TRBV:** first establish reference/numbering and real V calls.
   If paired data permit, fix donor-held-out partitions and receptor-only baseline
   before inspecting cell-state performance. No donor leakage or cell-level
   pseudoreplication. Insufficient data is an acceptable conclusion.

No formal RFU coloc without dense regional RFU association statistics and suitable
allele-aligned study LD. No calibrated Matos-only enrichment p-value without a
defensible RFU tested universe / exchangeability. No change to generic RFU assignment.

## Completed and reusable

- RfuWAS Data 1–4 preparation: 2,083 associations, 623 variants, 59 RFUs;
  623/623 REF checks passed previously and must not be repeated.
- Workbook SHA256: `42742f4a30548c1184f8de022a2adc7b6b9abd321e81a01ea740fad5527891d9`.
- eQTL archive: 1,743,867,610 bytes, verified MD5
  `e06e21a30576e6e271d974f922793ea6`; chr6/7 lookup and extraction manifests saved.
- Verified independent eQTL: 17 variants, 19 target/rank records, 15 genes, 39 RFUs.
- Matos workbook inside `PMC12870616_supplementary.zip`: verified MD5
  `2a8ed13f379946c1a7efce117c961201`. This is the XLSX hash, not its ZIP hash.
- Workbook Tables 1–6 respectively cover ChromBPNet scores, motifs, disrupted
  motifs, GWAS studies, QTL–GWAS colocalization, and ChromBPNet/GWAS CS summaries.
  Table 5 has 8,423 rows; no exact RFU-QTL variant intersections identified.
  Its GWAS PP.H4 values must not be relabelled molecular eQTL–caQTL PP.H4.
- Local RFU reference `~/ext/RFU-official/km5000noMax.Rdata` has 1,199,971
  CDR3-named cluster entries and 5,000 centers. No V-gene calls were found in that
  object; sequence motifs alone cannot establish exact V calls.
- Ensembl lookup identifies ENSG00000289938 and ENSG00000288882 as GRCh38 lncRNAs.
  They are non-TCR annotated targets, not established independent non-TCR mechanisms.

## Completed stages / final review

Stage 1 complete: caQTL **4,273,913,755 bytes**, MD5
`0a247f926008e7e7792a7689d327709a`. The full download completed without needing
assembly of the separately downloaded tail ranges. No download/reference check
needs repeating. All chr6/7 nominal and conditional members processed; partial
member SHA256 values match full extraction. New query family: 92,624 tests.
Complete results are in `results/final_pilot/`, preserving the partial checkpoint.
Conditional overlap: eQTL 17 variants / 19 target-rank records / 15 genes / 39 RFUs;
caQTL 32 variants / 41 target-rank records / 35 peaks / 39 RFUs; both at the same
variant 6 variants / 36 RFUs. Different-variant support within RFUs: 39 RFUs.

Persistent executable: `environment/bin/python`, resolving scRFU to
`/home/liuyuchen/Github/scRFU/src/scrfu/__init__.py`. Package freeze and R version
are in `logs/`; public interpretation source files are in `sources/interpretation_audit/`.

Stages 2/3 completed: workbook inspected and exact subsets cached; descriptive
Matos-only background computed without p-values; GSE190905 original held-out fits
completed. Their six-donor clone-weighted mean log loss increased from 1.852679
to 1.856996 with RFU. Original `rfu_state_validation_v1/` is preserved.
Candidate-RFU coverage is insufficient; study-to-study reference identity is
unverified. Full RFU locus statistics and source molecular-coloc outputs remain
the main formal-inference blockers.

Critical-review pass: the previous postprocessing process exited successfully,
including receptor joins and both figures. Its missing completion marker and
missing saved predictions were reproducibility gaps, not evidence of failed fits.
New runs mark `completion.json` running before work and complete only after all
required outputs have hashes. Consumers must check the manifest, not file existence.
The final dossier manifest and both audited prediction manifests were verified
against all listed outputs (9, 6 and 6 files). Recomputed saved-prediction scores
agree with fold metrics within 1e-12, and the original rerun matches exactly.
Ordinary donor generalization also completed: clone-weighted equal-donor mean
baseline/+RFU log loss 1.850996/1.855212, delta +0.004216; deterioration in all six
donors. No downloads or model fitting are outstanding.

## Claim-to-evidence ledger

### Frozen paired-data specification (before model fitting)

Use cached GSE190905 release: 6 donors, 27,655 TCR-bearing cells, 13 published
expression clusters. Final paper reports 7 donors / 57,738 cells; these are not
the dimensions of the cached release. Keep only correlation-threshold-qualified
RFUs (0.6), unambiguous primary TRBV/TRBJ and valid CDR3. No reassignment.
Require at least 4 donors with 100 eligible cells each. Leave one donor out;
exclude training V/J/CDR3 amino-acid identities observed in the held-out donor.
Retain evaluation classes with at least 20 training donor-clonotypes in at least
3 training donors; report excluded test coverage. RFU categories require at
least 10 training donor-clonotypes and 3 training donors; others pool as rare.
Baseline: one-hot TRBV, TRBJ, CDR3 length. Extended: baseline plus RFU category.
Fixed multinomial logistic regression, C=1, no hyperparameter search, no outcome-
guided feature selection. Training encoders and support rules use training only.
Primary weights: each donor V/J/CDR3 clone has total weight one, then balance
training donors. Sensitivity: cell weights, also balance training donors.
Evaluate weighted log loss and balanced accuracy separately in each held-out
donor; report paired donor differences and range, without calibrated p-values.
This is a single-cohort exploratory prediction comparison. Source clustering
used expression markers; exclusion of TCR genes from its HVG list is not documented.
Candidate RFU number compatibility with RfuWAS remains conditional until an
identical upstream reference checksum is available from that study.

### Review sensitivity specification (before fitting)

Preserve the original shared-receptor-excluded estimand. Refit it once in
`rfu_state_validation_audited_purged/` solely to save probabilities and verify
per-donor scores against the original; no model/settings changes. Also run
`rfu_state_validation_donor_generalization/` with shared receptors retained
across *different* donors. This addresses ordinary donor generalization and is
a separately reported sensitivity, not a replacement selected for its result.
All other filters, donor splits, regularization and weighting remain fixed.
Delta = extended log loss minus baseline; positive is worse. Audit paired cell
IDs/labels/weights and class-labelled probabilities. Six donors support descriptive
paired differences, not cell-level uncertainty or a population-wide negative claim.

| Claim | Evidence | Current status |
| --- | --- | --- |
| Exact RFU–eQTL association overlap | Prepared chr6/7 nominal and independent tables | Verified |
| Source conditional caQTL overlap | Verified chr6/7 members; 32 variants / 41 peak-rank records | Verified source selection, not LD-independent mechanisms |
| RFU molecular colocalization | Requires full RFU locus statistics | Not established |
| Matos gene–peak coloc at candidate | Requires released molecular signal-pair outputs | Not established |
| RFU 1415 specifically explained by TRBV28 | Competing TRBV18/TRBC1/TRBV7-7/lncRNA signals | Unsupported as an exclusive explanation |
| RFU adds information beyond V usage | Six-donor original and ordinary donor-generalization fits; low RFU coverage | No improvement demonstrated in either fixed design; not a universal negative |
| Data 4 independently validates RFU genetics | Same genetic-prediction framework | Not a valid claim |

## Checkpoint verification / next scientific step

```bash
# Read the explicit completion markers; do not rerun completed models or downloads.
cat ~/data/scRFU_regulatory_20260909/results/final_pilot/completion.json
cat ~/data/scRFU_regulatory_20260909/results/rfu_state_validation_audited_purged/completion.json
cat ~/data/scRFU_regulatory_20260909/results/rfu_state_validation_donor_generalization/completion.json
```

Next scientific dependency: dense regional RFU-QTL summary statistics with
verified effect alleles and appropriate study LD. See the unsent data-request
draft in the pilot report. No formal RFU colocalization is currently ready.

## Final checkpoint validation

- Focused regulatory/preparation/reference/audit tests: **104 passed**.
- Focused prediction/completion tests after the title-only fix: **9 passed**.
- Full suite after final code changes: **466 passed, 4 skipped, 32 warnings**
  (22.09 s). Command: `NUMBA_CACHE_DIR=/tmp/scrfu-regulatory-numba-cache
  environment/bin/python -m pytest -q -o addopts=''` from the repository, with
  `environment/bin/python` resolved beneath the external workspace.
- `ruff check .`: passed. `ruff format --check .`: **218 files already formatted**.
  `git diff --check`: passed. The two formatted documentation examples have
  identical Python ASTs before/after formatting.
- Final evidence/donor manifests: **21 output hashes verified** across three
  completed stages. Saved prediction audits and figure labels agree with the
  report; the wrapped title is legible. Core `src/scrfu/regulatory.py` unchanged.
- Logs are retained in the external workspace `logs/validation_checkpoint/`.
  This checkpoint is complete; no further analysis should run automatically.
