# Submission gate

Freeze date: 2026-08-24. The machine-readable checklist is
[`submission_gate.json`](submission_gate.json). “Complete” requires executed,
audited evidence or a repository-verifiable artifact; a software function is
not evidence. Runtime source tables remain external by design.

| Gate | Assessment | Evidence | Blocker / exact next action |
|---|---|---|---|
| Original-RFU parity | Complete | adversarial real official-reference parity and immutable hashes | Preserve artifacts and score tolerance. |
| Novelty beyond wrapping | Strong partial | deterministic scaling, robustness, transfer, held-out evaluation, phenotype linkage, and comparators | Run governed longitudinal validation or narrow the longitudinal manuscript claim. |
| Governed longitudinal validation | Blocked | frozen methods and exact runtime command only | Configure the three explicit paths; do not search for or move private inputs. |
| VDJdb antigen evidence | Blocked | release-pinned acquisition and frozen null design only | Configure the 2026-06-03 reference and run, or remove the antigen claim and use Figure 4 fallback. |
| Two independent public TCR datasets | Complete for technical transfer | Wells, GSE190905, and GSE157007 with distinct roles | An additional aging development cohort is optional unless an aging claim is added. |
| Completely held-out dataset | Complete | preregistered GSE157007 evaluation | Retain cross-sectional one-sample-per-donor caveat. |
| BCR | Deliberately removed from first manuscript | scope freeze | Do not implement or present BCR. |
| Meaningful comparators | Complete for executed tasks | RFU and conventional sample-level representations with identical samples/candidates | Report complementarity; edit-distance remains unavailable above the frozen 2,000-sequence limit. |
| Runtime and memory | Complete for bounded scale | Wells 1k/10k/25k fresh/resume/serial/parallel records | Report hardware and do not extrapolate to all atlas cells. |
| Methods and claim freeze | Complete | methods freeze, Results claims, figure plan, and claim audit | Apply documented change control only. |
| Installation/API documentation | Strong partial | API freeze, README, examples, CI wheel smoke | Add and clean-wheel-test an end-to-end tutorial and final public API presentation. |
| Public test dataset | Weak partial | deterministic synthetic test fixtures | Designate a licensed small public or generated tutorial fixture with expected hashes. |
| CI | Complete | Python 3.10–3.12 tests/Ruff/build/wheel smoke | Require green CI on the exact release commit. |
| Open-source license | Complete | MIT license | Verify inclusion in release artifacts. |
| Citation metadata | Blocked | no `CITATION.cff` | Add after 0.4.0 approval without inventing DOI or ORCID. |
| Versioned release and DOI | Blocked | release gap audit | Reconcile 0.1.0 source with historical v0.2.0/v0.3.0 tags as 0.4.0, then archive only after approval. |
| Reviewer-shareable code/data | Strong partial | public manifests, audit, inventory, and acquisition instructions | Prepare the exact release and disclosure-reviewed source-data deposit. |
| Complete source tables | Strong partial | Figures 1 and 3 sources audited | Figure 2 is blocked; primary Figure 4 is blocked unless its fallback is selected. |
| Manuscript skeleton | Strong partial | demonstrated prose and explicit blocked markers | Resolve or remove blocked Results claims, add references, and adapt to journal format. |
| Reproducibility report | Strong partial | claim audit, Month 2 report, methods freeze, and source inventory | Append governed/optional antigen evidence and exact release provenance if generated. |

## Minimum blockers for Cell Reports Methods

1. **Scientific scope:** the manuscript's longitudinal claim requires the real
   governed repeated-measures run. If that input cannot be supplied, the title,
   abstract, claim structure, and Figure 2 must be narrowed rather than filled
   with simulated evidence.
2. **Figure 4 decision:** VDJdb is required only if antigen-evidence coherence is
   retained. It is not a blocker to the technical transfer story if the claim is
   removed and the prespecified robustness/transfer fallback is used.
3. **Submission artifacts:** final disclosure-reviewed source tables, an
   end-to-end tutorial/public fixture decision, `CITATION.cff`, a clean 0.4.0
   release-candidate install test, and the final manuscript/references are still
   required before submission.
4. **Release/archive:** a versioned release and DOI are required for the final
   public submission package, but they are downstream of scientific and release
   approval and are not authorized by this gate.

Optional additional cohorts, speculative methods, and BCR are not minimum
blockers for the current central methods story. No gate here authorizes
publication, tagging, or private-data movement.
