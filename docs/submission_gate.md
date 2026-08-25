# Submission gate

Freeze date: 2026-08-24. The machine-readable checklist is
[`submission_gate.json`](submission_gate.json). “Complete” requires executed,
audited evidence or a repository-verifiable artifact; a software function is
not evidence. Runtime source tables remain external by design.

| Gate | Assessment | Evidence | Blocker / exact next action |
|---|---|---|---|
| Original-RFU parity | Complete | adversarial real official-reference parity and immutable hashes | Preserve artifacts and score tolerance. |
| Novelty beyond wrapping | Complete for the independent methods/software scope | deterministic full-atlas scaling, robustness, transfer, held-out evaluation, public paired structure, phenotype linkage, VDJdb evidence, and comparators | Retain the limited technical/transfer interpretation. |
| Public longitudinal methods validation | Complete for paired donor structure | GSE190905 public pre/post within/between and leave-one-timepoint-out retrieval | Temporal dynamics/compartment claims require a deeper public cohort and are not part of the current gate. |
| Private longitudinal cohort | Deliberately removed as a dependency | independent tool evidence is public-only | Do not use or search for private cohort data for this paper/release. |
| VDJdb antigen evidence | Complete | pinned 2026-06-03 reference; 24 sensitivity cells and four 1,000-permutation null designs per supported cell | Report annotation coherence and sparse strict-match cases; never claim antigen specificity. |
| Two independent public TCR datasets | Complete for technical transfer | Wells, GSE190905, and GSE157007 with distinct roles | An additional aging development cohort is optional unless an aging claim is added. |
| Completely held-out dataset | Complete | preregistered GSE157007 evaluation | Retain cross-sectional one-sample-per-donor caveat. |
| BCR | Deliberately removed from first manuscript | scope freeze | Do not implement or present BCR. |
| Meaningful comparators | Complete for executed tasks | RFU and conventional sample-level representations with identical samples/candidates | Report complementarity; edit-distance remains unavailable above the frozen 2,000-sequence limit. |
| Runtime and memory | Complete | Wells 1k/10k/25k benchmarks plus full 610k-cell receptor-only run and downstream stress passes | Report hardware and distinguish assignment from downstream timing. |
| Methods and claim freeze | Complete | methods freeze, Results claims, figure plan, and claim audit | Apply documented change control only. |
| Installation/API documentation | Complete for release candidate | API freeze, README, runnable tutorial, examples, CI wheel/tutorial smoke | Hosted API presentation is a strong recommendation, not a blocker. |
| Public test dataset | Complete | packaged fully synthetic tutorial fixtures with expected hashes | Preserve fixture hashes and synthetic provenance. |
| CI | Complete | Python 3.10–3.12 tests/Ruff/build/wheel smoke | Require green CI on the exact release commit. |
| Open-source license | Complete | MIT license | Verify inclusion in release artifacts. |
| Citation metadata | Complete for candidate | `CITATION.cff` contains author, repository, license and 0.4.0rc1 without DOI/ORCID | Change version only with final release approval. |
| Versioned release and DOI | Blocked by approval | coherent 0.4.0rc1 candidate; historical tags untouched | Review, approve, then finalize 0.4.0 and archive in a later authorized session. |
| Reviewer-shareable code/data | Strong partial | public manifests, audit, inventory, and acquisition instructions | Prepare the exact release and disclosure-reviewed source-data deposit. |
| Complete technical source tables | Strong partial | full Wells and VDJdb source summaries are external; public transfer/held-out sources audited | Assemble the final public disclosure-reviewed deposit; deeper longitudinal dynamics are outside current scope. |
| Manuscript skeleton | Strong partial | demonstrated prose and explicit blocked markers | Resolve or remove blocked Results claims, add references, and adapt to journal format. |
| Reproducibility report | Strong partial | claim audit, Month 2 report, methods freeze, and source inventory | Append governed/optional antigen evidence and exact release provenance if generated. |

## Minimum blockers for Cell Reports Methods

1. **Submission artifacts:** final disclosure-reviewed public source tables,
   clean final 0.4.0 installation validation, and the final
   manuscript/references are still required before submission.
2. **Scientific wording:** claims must remain technical/transfer-focused;
   GSE190905 is a six-donor paired methods demonstration and VDJdb matches are
   annotation coherence, not specificity.
3. **Release/archive:** a versioned release and DOI are required for the final
   public submission package, but they are downstream of scientific and release
   approval and are not authorized by this gate.

Optional additional cohorts, speculative methods, and BCR are not minimum
blockers for the current central methods story. No gate here authorizes
publication, tagging, or private-data movement.
