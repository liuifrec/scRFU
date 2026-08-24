# Manuscript claim audit

Frozen audit date: 2026-08-24. Runtime paths are relative to
`SCRFU_MONTH2_OUTDIR`. The audited runtime root used for this review was the
external Month 2 output directory; no runtime table is copied into Git.

## Shared provenance

| ID | Provenance |
|---|---|
| R1 | Official RFU commit `ff7a1ea6aca28444a7ff8ec7ef09745e7341f99d`; `RFU.R` SHA256 `92c2faa33f2e7f60d6470ad7dfd653eb0ca54859aa198ed53ad0c876aae8640b`; trimer SHA256 `820fb71428913974e994543cbbbaa591f54355cdb73d651f068f0e691a6cdffd`; centroid SHA256 `64783074360edeca84b3f49ab9d682263e954752fc09bb208edaabde0fd82553`; standard mode; threshold 0.6. |
| P1 | Parity original-input SHA256 `9da924a4c003b74aa48e4e8fc21260f72e37b1e197d9521da015ea2be39077f3`; canonical scRFU-input SHA256 `5a4f781ee67ff9d680e0970decb416b6ed7741da9ca4bd724849efb10de5e81b`; seed 20260824. |
| W1 | Wells source fingerprint `212852f0d80c79184f7625e379c15d66f41c93d7a65c3b6cde82c2ecde371c8a` using `sha256-size-first-middle-last-1MiB-v1`; selected-cell hashes: 1k `ce40e392f142c19e75d9e7160bb930b0653e76536d83ae6597b8807c8ca64a68`, 10k `1433550c5dda719f6ea83eb49345083ceef817ccf87dfc4dd66ce442fd9af658`, 25k `279109ec6dc4b019a964cf8e56abb367a4c40f6856bb50f9e1f6de36e9780eda`; scaling-input hashes 1k `27b962d0bc925c41392507269d4e242a81d364116141b284301e7f7e034ad2f5`, 10k `cd75d031697ec90adc8759531a2cbb7a329c590177e740c34d74b616172df798`, 25k `e10eaf772ff1a98a595517b9677fd397fce066a20d3fc8c483348c4a8e01ae06`; seed 20260824. |
| G1 | GSE190905 receptor SHA256 `85a9db858306b68d328b82b53b99ed38f3f4228d18c06a446c4d461c813e47c1`; metadata SHA256 `568297bbba944323b860cc3df89e7d1377c04bcab93246eca94749886c617537`; prepared H5AD SHA256 `c3d966f41d93962aea0124fc2d7022c4e26b95d8e8554e69721e9a952be492b3`. |
| H1 | GSE157007 canonical 17-file receptor-manifest SHA256 `6e311d8c537da80f1aa023e06c5708c58ebbf4e4ce902415d98074db1ca75732`; family-metadata SHA256 `61056e7766a7cca0bbe8d24b99e09a31cbc9588519b488332a45e525c75fd95f`; preregistration fixes productive primary TRB, cell weighting, standard mode, threshold 0.6, and nearest/threshold-qualified policies. |

All executed analyses below used scRFU 0.1.0 on 2026-08-24 and R1 unless a
row says otherwise. “Technical” covers software equivalence, runtime, memory,
coverage, and perturbation behavior. “Descriptive biological” covers observed
sample/donor/phenotype structure without population inference. No audited claim
is inferential.

## Quantitative statement ledger

| Manuscript-facing quantitative statement | Runtime source and SHA256 | Input / parameters | Biological unit | Class | Audit decision |
|---|---|---|---|---|---|
| Parity fixture: 12 total rows, 10 eligible, 2 ineligible, 5 unique eligible CDR3s; 10/10 RFU IDs matched; label, threshold, reconstruction, and ordering mismatches all 0; maximum score error 0 at tolerance `1e-12`. | `original_rfu_parity/summary.json`, `4adc57d3415e6377cee4dcfe2be450a906c688fbd1b80e9a40dff76f0ff943e5` | P1; R1; threshold 0.6 | receptor row and unique CDR3 | Technical | Verified. |
| Parity execution took 3.924 s; official child peak RSS was 352,116 KB. | same parity summary | P1; R1 | one bounded execution | Technical | Verified; hardware-specific. |
| Wells source contained 610,429 cells, 24 donors, 17 cell types, 10 tissues, and 247 libraries. | `wells_bounded/wells_sampling_summary.tsv`, `2b2898788c782d0cd15c3bd17481222a8c86b9df1a2d89abca0deac6cab8c197` | W1; selective observation/receptor reads | source cell / donor / annotation category | Descriptive biological | Verified; no expression matrix was read. |
| Random 1k/10k/25k subsets contained 499/4,861/12,415 productive TRB rows and 453/4,130/10,038 unique queries. | sampling table plus `wells_bounded/figure1_scaling_source.tsv`, `100bcb6c8dff56256e7e91df3cf7a5598cf30e44717a7c3769e1bcf4e7c276e4` | W1; deterministic sampling; source order retained | source cell / receptor row / unique CDR3 | Technical | Verified. |
| Donor×tissue×cell-type stratum total variation was 0.0772 for random 25k and 0.0115 for proportional stratified 25k. | Wells sampling table | W1; hash-ranked proportional allocation | sampled cell distribution | Technical | Verified; neither subset is called representative. |
| Unchunked serial RFU times were 3.920/18.141/40.276 s at 1k/10k/25k. | Wells scaling table | W1; unchunked, one process worker | one RFU execution | Technical | Verified; hardware-specific. |
| At 25k with 5,000-query chunks, one worker took 45.650 s and two workers 25.859 s. | Wells scaling table | W1; three chunks; process executor | one RFU execution | Technical | Verified; primary scaling contrast. |
| At 25k with 1,000-query chunks, one worker took 63.809 s and two workers 36.102 s. | Wells scaling table | W1; 11 chunks; process executor | one RFU execution | Technical | Verified; secondary scaling contrast. |
| Resumed 25k runs took 0.209–0.222 s and reused every completed chunk. | Wells scaling table | W1; 1,000/5,000-query chunks; one/two workers | one resumed execution | Technical | Verified. |
| Wells 25k Python peak RSS was 482,296 KB; maximum reported backend RSS was 1,802,608 KB. | Wells scaling table | W1 | one execution process / child backend | Technical | Verified; peak values are hardware/runtime-specific. |
| Unique-query/input ratios were 0.908/0.850/0.809 for Wells 1k/10k/25k. | Wells scaling table | W1; exact CDR3 deduplication | receptor row / unique CDR3 | Technical | Verified. |
| Wells 25k observed 3,963 nearest RFUs and 3,381 threshold-qualified RFUs; 9,552/12,415 cells passed threshold (76.94%). | `wells_bounded/single_cell_analysis/assignment_policy_summary.tsv`, `f3de40b04831e523b26c1c26700c49a7f9707383464c2125120048321d2ad120`; scaling table | W1; nearest and threshold-qualified policies | receptor-bearing cell / RFU | Technical | Verified. |
| All meaningful serial/parallel/chunk/resume/shuffle runs share assignment SHA256 `3c57b4e983fa1dc9150d99ff3dfa6f6b6cdf7c5db86675cc93e287540f5317b3` at 25k. | Wells scaling table | W1; score tolerance `1e-12` | receptor row | Technical | Verified. |
| Nearest-policy cell-subsampling median cosine was 0.526/0.730/0.885/1.000 and median Spearman 0.510/0.713/0.873/1.000 at 25/50/75/100%. | `wells_bounded/figure1_robustness_source.tsv`, `74ea9f21c5c9005a070b9e75b476558aeaf9c342ba09a114d2be8f05e435ceb1` | W1; three fixed seeds; full bounded reference | RFU abundance vector | Technical | Verified and corrected from prior prose transcription. |
| Nearest-policy phenotype-coupling cosine means were 0.506/0.708/0.863/1.000 and dominant-phenotype agreement 0.637/0.735/0.853/1.000 at 25/50/75/100% cell subsampling. | `wells_bounded/phenotype_coupling_stability.tsv`, `205a4da7606517ad571f7a0ce0f6a356bef42ca480f6a3a87b772742d0c5d385` | W1; three fixed seeds; cell type phenotype | RFU×phenotype coupling profile | Descriptive biological | Verified; no cell-level inferential test. |
| Wells nearest/threshold analyses covered 16 cell types, 21 receptor-bearing donors, and 209 libraries. | Wells assignment-policy summary | W1 | cell type / donor / library | Descriptive biological | Verified; counts differ from whole-atlas metadata because only receptor-bearing rows enter. |
| GSE190905 contained 27,655 receptor rows, 16,391 unique CDR3s, 6 donors, and 12 paired time samples. | `public_data/GSE190905/validation/run_manifest.json`, `00eb04ebe972d2425187520a896ae654ddc9781e7d4210523720043c6a8c5378` | G1; primary TRB; 5,000-query chunks; two workers | receptor cell / patient-time sample | Descriptive biological | Verified. |
| GSE190905 threshold coverage was 21,657/27,655 (78.31%); 4,465 nearest and 3,940 threshold-qualified RFUs were observed. | same GSE190905 manifest | G1; R1 | receptor cell / RFU | Technical | Verified. |
| GSE190905 RFU execution took 45.993 s; backend peak RSS was 946,580 KB. | same GSE190905 manifest | G1; 5,000-query chunks; two workers | one execution | Technical | Verified; hardware-specific. |
| In 12 leave-one-timepoint-out queries, nearest-RFU cosine retrieval was top-1 0.667, top-3 0.833, MRR 0.764; exact-CDR3 cosine was 0.667/0.750/0.771 and clonotype cosine 0.667/0.833/0.756. | `cross_cohort/gse190905_donor_retrieval_comparators.tsv`, `d2f79d96b4e0ea39bef170d1fe2a37ff9b5aec95e314f70c993f2d23f6ab468f` | G1; identical candidate sets; leave one timepoint out | patient-time query | Descriptive biological | Verified; six-patient methods demonstration, no population inference. |
| GSE190905 nearest-RFU mean cosine was 0.491 within donor versus 0.063 between donor; threshold-qualified RFU was 0.500 versus 0.055. | `cross_cohort/gse190905_within_between_comparators.tsv`, `8911ba5abd980c321d17b13dbb5f1cc6028748e5d8c08718131c74d13021b8e5` | G1; 6 within-donor and 60 between-donor pairs | patient-time sample pair | Descriptive biological | Verified; no donor-aware inferential test was run. |
| Held-out GSE157007 contained 60,125 productive primary-TRB receptor rows, 43,082 unique CDR3s, 17 donors, and 4 age groups. | `public_data/GSE157007/heldout_validation/run_manifest.json`, `c2180b82ec514aa311a9fda5667ca2f900dbb7be34b0dd463d6d5d4f036c6461` | H1; R1; 5,000-query chunks; two workers | receptor cell / donor | Descriptive biological | Verified against preregistration. |
| Held-out threshold coverage was 76.88%; 4,898 nearest and 4,618 threshold-qualified RFUs were observed. | `public_data/GSE157007/heldout_validation/transfer_summary.tsv`, `5a02dbcb956de7047b1bdf83dcfcf912887ee817b11def9bb88623a3f0525a56` | H1; R1 | receptor cell / RFU | Technical | Verified. |
| Held-out nearest-RFU mean sample-vector cosine was 0.936 at 50% and 0.976 at 75% deterministic subsampling; threshold-qualified values were 0.922 and 0.970. | `public_data/GSE157007/heldout_validation/subsampling_stability_summary.tsv`, `b5cbcbb1edc3ee6fb1806092833a66af695bb009391f21f8d87d8d5964f7c866` | H1; three seeds; 51 sample-seed observations per fraction | biological-sample vector | Technical | Verified; one sample per donor prevents donor retrieval. |
| Conventional comparator dimensions were: GSE190905 exact CDR3 16,391, clonotype 10,407, V 49, J 13, length 15; GSE157007 exact CDR3 43,082, namespaced clonotype 45,899, V 49, J 13, length 16. | GSE190905 comparator status, `16eb1beb02c68d3d75c1236140ac5d630fd1c3cefb33aee813b54793e0dfc36e`; GSE157007 comparator status, `813b83df82a00d214fff42580b26ce6dc70ffd0eb10cf495a965ac25760c73a3` | G1/H1; sample-level count matrices | biological sample / comparator feature | Technical | Verified; edit-distance skipped above frozen 2,000-sequence limit. |

## Public cohort-design statements

These are public-source metadata rather than runtime results. PRJNA602091's 30
adults and two visits averaging 9.2 years, GSE190905's six-patient paired design,
GSE157007's 17-donor age/frailty design, and GSE158848's 86 repertoires from five
donors remain traceable to the publication/GEO links in
[`public_dataset_candidates.md`](public_dataset_candidates.md). They are cohort
descriptions, not scRFU outcome claims.

## Excluded or blocked quantitative statements

- GSE157007 fresh end-to-end time was observed but not retained in a dedicated
  immutable summary before the resume manifest replaced the top-level timing.
  Chunk manifests remain, but that end-to-end number is excluded from figures
  and Results prose.
- VDJdb match, coherence, ambiguity, and permutation statistics do not exist
  because no release-pinned runtime reference is configured.
- Governed longitudinal cohort counts, similarities, retrieval, dynamics,
  compartment contrasts, and resampling statistics do not exist because the
  three required runtime paths are unset.
- No cell-level or sequence-level p-value is treated as biological inference.
