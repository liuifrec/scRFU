# scRFU figure prototype QA

This audit covers the 22 main panels generated from frozen evidence. “Standalone”
means that the plotted object, analysis unit and caveat are legible without relying
on Results prose. It does not mean that a full figure legend is unnecessary.

## Figure 1

| Panel | Exact supported claim | Source support | Analysis unit obvious? | Comparator fair? | Overclaim risk and control | Standalone? |
|---|---|---|---|---|---|---|
| 1A | A frozen reference maps exact receptor queries into a shared RFU vocabulary while exact identity remains available. | Yes; conceptual description only. | Yes: receptor sequence/feature. | Yes; exact identity is retained, not presented as inferior. | Medium; shared vocabulary is not biological equivalence. The panel says “representation,” not “functionally identical.” | Yes. |
| 1B | Chain-aligned output preserves all AIRR chain positions and limits assignment to eligible TRB chains. | Yes; exact native-schema tests. | Yes: AIRR chain. | Not applicable. | Low; explicitly labelled as a software invariant. | Yes. |
| 1C | RFU records and portable provenance coexist with AIRR data without expression access. | Yes; storage schema and X-independent tests. | Yes: observation and nested chain. | Not applicable. | Low; no backed-mode claim. | Yes. |
| 1D | Cell summaries are opt-in and ambiguity-aware under conflicting eligible chains. | Yes; exact policy tests. | Yes: cell after chain-level assignment. | Yes; agreement and conflict are both shown. | Low; no biological preference among policies is implied. | Yes. |
| 1E | The public Scirpy smoke produced exact native/table parity for 2,931 eligible TRB chains. | Yes; public run summary and manifest. | Yes: eligible TRB chain. | Yes; same backend/input sequences. | Low; labelled as a 3,000-cell interoperability smoke and Linux RFU run. | Yes. |

## Figure 2

| Panel | Exact supported claim | Source support | Analysis unit obvious? | Comparator fair? | Overclaim risk and control | Standalone? |
|---|---|---|---|---|---|---|
| 2A | The bounded adversarial fixture exactly matches official RFU IDs, labels, scores, threshold state and row reconstruction. | Yes; audited parity ledger and hashed external table. | Yes: eligible receptor row. | Yes; official `AssignRFUs()` is the reference. | Medium if generalized to all possible inputs; the bounded-fixture caveat remains explicit. | Yes. |
| 2B | Targeted receptor-only execution reduces 303,088 productive rows to 192,675 unique queries and processes the full Wells atlas in 8.5 minutes at 3.42 GB peak RSS. | Yes; sealed full-run evidence. | Yes: cells, rows and unique queries are separately named. | Yes; the funnel compares stages, not methods. | Medium; performance is host-specific and not expression-inclusive. | Yes. |
| 2C | Native execution scales from 25k to 250k cells with zero assignment mismatches against the table path. | Yes; three native source tables. | Yes: selected cell/eligible chain. | Yes; identical sequences and frozen backend. | Low; only one Linux host and bounded native objects are claimed. | Yes. |
| 2D | Compatible caches markedly shorten repeat execution while tested serial/parallel/chunk/order configurations preserve assignments. | Yes; fresh/resume manifests and audited invariant ledger. | Yes: execution configuration. | Yes; fresh and cached runs share inputs/reference/parameters. | Medium; incompatible caches are not claimed reusable. | Yes. |

## Figure 3

| Panel | Exact supported claim | Source support | Analysis unit obvious? | Comparator fair? | Overclaim risk and control | Standalone? |
|---|---|---|---|---|---|---|
| 3A | The observed RFU vocabulary is smaller than the exact-CDR3 vocabulary in all four public datasets. | Yes; frozen compression table. | Yes: dataset feature vocabulary. | Yes; counts use the same eligible receptors per dataset. | High if equated with better biology; title and caveat describe compression only. | Yes. |
| 3B | RFUs aggregate dataset-dependent numbers of distinct CDR3 sequences. | Yes; per-RFU frozen counts. | Yes: RFU. | Yes; curves use identical nearest-assignment semantics. | High if interpreted as proven functional equivalence; explicitly disclaimed. | Yes. |
| 3C | RFU occupies an intermediate dimension/sparsity regime relative to exact identity and coarse summaries in GSE190905. | Yes; genuine Scirpy comparator table. | Yes: sample-by-feature matrix. | Yes; identical 12 samples. | Medium; neither low dimension nor low sparsity is labelled superior. | Yes. |
| 3D | One frozen threshold/reference yields approximately 77% coverage across four public datasets. | Yes; frozen compression table. | Yes: eligible sequence/chain. | Yes; same threshold/reference. | High if treated as calibrated OOD probability; explicitly disclaimed. | Yes. |
| 3E | RFU identities are shared across datasets much more frequently than exact CDR3 identities. | Yes; exact set-overlap table. | Yes: dataset pair and feature identity. | Yes; exact-set Jaccard is used for both representations. | High; shared feature space is not biological equivalence and the annotation says “represented/shared.” | Yes. |

## Figure 4

| Panel | Exact supported claim | Source support | Analysis unit obvious? | Comparator fair? | Overclaim risk and control | Standalone? |
|---|---|---|---|---|---|---|
| 4A | In GSE190905, same-donor two-visit samples are more similar than between-donor samples for RFU, exact CDR3, Scirpy clonotype and V/J representations. | Yes; frozen pairwise summary. | Yes: patient-time sample pair. | Yes; fixed samples and cosine metric. | Medium; pairwise means are descriptive and the cohort has only two visits. | Yes. |
| 4B | RFU, exact CDR3, Scirpy clonotype and V/J retrieve donors comparably on some endpoints, while exact/VJ lead RFU on MRR and coarse summaries differ. | Yes; frozen retrieval summary. | Yes: held-out patient-time sample. | Yes; identical candidate sets. | Low; all methods and non-RFU-leading results remain visible. | Yes. |
| 4C | Stability under 50%/75% subsampling varies by representation, with coarse summaries often most stable. | Yes; fixed-seed frozen summary. | Yes: sample-seed representation. | Yes; same fractions/seeds. | Medium; stability can reflect information loss, which is stated. | Yes. |
| 4D | The preregistered held-out cohort maps at 76.88% threshold coverage and retains high RFU-vector cosine under 50%/75% subsampling. | Yes; sealed held-out manifest and audited source-table hashes. | Yes: receptor and sample vector. | Yes; nearest and threshold policies were preregistered. | Medium; one sample per donor means no held-out donor-retrieval claim. | Yes. |

## Figure 5

| Panel | Exact supported claim | Source support | Analysis unit obvious? | Comparator fair? | Overclaim risk and control | Standalone? |
|---|---|---|---|---|---|---|
| 5A | Abundant RFUs retain heterogeneous single-cell phenotype profiles under nearest and threshold-qualified policies. | Yes; full Wells coupling tables. | Yes: RFU × cell type. | Yes; policies share the same abundance-selected RFUs. | High if RFUs are outcome-selected; selection is explicitly total-abundance only. | Yes. |
| 5B | Bounded Wells phenotype-coupling profiles become increasingly concordant with the full reference as retained-cell fraction increases. | Yes; audited frozen stability table/hash. | Yes: RFU × phenotype profile. | Yes; fixed seeds and fractions. | Medium; descriptive stability only, without cell-level inferential tests. | Yes. |
| 5C | Native chain-aligned output can link CDR3-based RFU identity to exact CDR3 and strict CDR3+V VDJdb evidence without chain expansion. | Yes; native linkage tables and invariant manifest. | Yes: AIRR chain. | Yes; both match keys are shown separately. | High if identities are conflated; the panel explicitly keeps RFU identity exact-CDR3 based. | Yes. |
| 5D | Distinct Wells sequences grouped by RFU show higher external same-antigen annotation coherence than four prespecified size-preserving nulls. | Yes; 1,000-permutation frozen null table. | Yes: distinct matched RFU CDR3 sequence. | Yes; four named null controls are displayed. | High; this is external annotation coherence, not antigen specificity or prediction. | Yes. |

## QA disposition

- No main panel requires language stronger than its frozen evidence.
- Figure 3 is the strongest conceptual centerpiece but needs the representation-property
  caveats in its legend.
- Figure 5 is scientifically bounded by sparse VDJdb exact-match coverage and descriptive
  phenotype coupling; both limitations are visible rather than hidden.
- Figure 1 is intentionally schematic-heavy. Its empirical anchor is the public
  `wu2020_3k` parity panel.
- No new experiment is required for these figures or the frozen central claim.
  GSE345124 remains optional unless the scope expands to deep longitudinal dynamics.
