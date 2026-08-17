# Submission gate

The machine-readable checklist is [`submission_gate.json`](submission_gate.json).
Only repository evidence can support a `complete` status. Real-data conclusions,
private inputs, and planned analyses are not evidence by themselves.

| Gate | Status | Evidence | Main blocker / next action |
|---|---|---|---|
| Original-RFU parity | Partial | integration tests | Produce reviewed frozen-artifact parity tables. |
| Novelty beyond wrapping | Partial | methodological novelty and implemented methods | Validate longitudinal and transfer claims independently. |
| Longitudinal compartment results | Blocked | synthetic longitudinal module only | Run prespecified governed analysis outside this repository. |
| Two independent TCR datasets | Blocked | dataset plan | Verify and execute suitable public cohorts. |
| Completely held-out dataset | Blocked | held-out plan and manifest API | Register before access/evaluation. |
| BCR substantive or removed | Complete for TCR-only scope | BCR design gate | Keep omitted unless every gate passes. |
| Meaningful comparators | Partial | comparator interface | Run frozen comparisons on real cohorts. |
| Runtime and memory | Partial | benchmark plan/manifests | Run bounded scaling series. |
| Installation/API docs | Partial | API freeze and README | Complete tutorials and release-version review. |
| Public test dataset | Partial | offline synthetic examples | Decide on redistributable small public fixture. |
| CI | Complete | workflow | Monitor all supported Python jobs. |
| Open-source license | Complete | MIT license | Retain in all artifacts. |
| Versioned release and DOI | Blocked | changelog only | Approve release candidate, then archive. |
| Reviewer-shareable code/data | Partial | public-input workflow design | Assemble public manifests and instructions. |
| Complete source tables | Blocked | figure plan | Generate source data for every panel. |
| Reproducibility report | Blocked | this checklist | Compile only after final workflows and artifacts pass. |

No gate here authorizes publication, tagging, or private-data movement.
