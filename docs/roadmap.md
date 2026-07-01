# Roadmap

## v0.1: RFU Backend and AnnData Attachment

- Integrate with an external upstream RFU repository without vendoring RFU code
  or data.
- Extract AnnData/scirpy-style AIRR features needed for TRB RFU assignment.
- Attach RFU labels, scores, and provenance into AnnData.

## v0.2: Summary, Aggregation, Plotting, and Export

- Summarize per-cell RFU assignment rates and score distributions.
- Aggregate RFU labels by user-provided metadata groups.
- Provide lightweight matplotlib plotting utilities.
- Export group-by-RFU matrices for downstream analysis.

## v0.3: Synthetic Demo

- Provide a small synthetic scirpy-style AnnData example.
- Exercise extraction, annotation, summary, aggregation, plotting, and export
  without requiring R, RFU files, internet access, or real datasets.

## v0.4 Planned: Real-Data Demonstration

- Add reproducible workflow templates for user-provided public single-cell
  RNA/TCR datasets.
- Demonstrate assignment summaries and RFU abundance matrices on at least one
  real dataset.
- Keep data download and storage outside the repository.

## v0.5 Planned: Manuscript Workflow and Reproducible Figure Scripts

- Add scripts for producing manuscript-oriented summary tables and figures from
  user-provided inputs.
- Document expected input AnnData fields and output artifacts.
- Keep generated figures and data-heavy outputs out of version control unless
  they are small, synthetic, and CI-safe.

## Future: B-Cell/ALFU Experimental Backend

- Explore an experimental antibody-lineage functional unit design for BCR data.
- Treat this as a future extension rather than a current scRFU capability.

## Future: Radiation-Specific Workflows

- Provide templates for radiation-associated PBMC datasets where appropriate
  metadata are available.
- Avoid biological claims until validated on real datasets with appropriate
  controls.

## Future: ScGeo/Repertoire-State Topology Integration

- Consider future integration with repertoire-state topology methods after the
  RFU single-cell workflow is stable and documented.
