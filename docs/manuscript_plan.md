# Manuscript Plan

## 1. Working Title Options

- scRFU: single-cell integration of repertoire functional unit analysis in AnnData workflows
- scRFU enables repertoire functional unit annotation of single-cell TCR data
- Functional-unit annotation of adaptive immune repertoires in single-cell immune profiling

## 2. Core Contribution

scRFU bridges an established upstream RFU repertoire functional unit method with
single-cell RNA and TCR workflows. It connects 10x/scRNA+TCR data structures,
AnnData/scirpy-style AIRR storage, per-cell RFU annotation, group-level RFU
aggregation, visualization, export, and provenance tracking.

The main contribution is an integration framework: scRFU makes sequence-derived
RFU assignments usable alongside transcriptomic cell states, sample metadata,
cell annotations, and downstream scverse analysis objects.

## 3. What scRFU Is

- An integration framework for applying RFU annotation in single-cell immune
  profiling workflows.
- A Python wrapper around an external upstream RFU backend.
- An AnnData-compatible annotation layer for storing RFU labels and scores in
  `adata.obs`.
- A lightweight analysis layer for RFU summaries, group-level aggregation,
  plotting, export, and provenance.

## 4. What scRFU Is Not

- scRFU is not a new RFU algorithm.
- scRFU is not a replacement for the original RFU method or implementation.
- scRFU is not yet a full BCR/ALFU framework.
- scRFU is not yet a radiation biology discovery paper unless real-data
  validation is completed.

## 5. Proposed Manuscript Figures

### Figure 1: Workflow Schematic

10x VDJ/scirpy AnnData to TRB CDR3aa/TRBV extraction, upstream RFU execution,
`adata.obs` annotation, aggregation, plotting, export, and provenance.

### Figure 2: Synthetic/Demo Workflow Output

- RFU assignment summary.
- Group-by-RFU heatmap.
- RFU score histogram.

### Figure 3: Real Public Dataset Demonstration

- Assignment rate by sample or cell type.
- RFU abundance heatmap.
- UMAP or metadata-linked RFU distribution if available.

### Figure 4 Optional: Radiation PBMC Demonstration

- RFU composition differences between condition, time, or dose groups.
- Enrichment or deviation summary.

This figure should be included only if the data pipeline succeeds and the
analysis supports the displayed results.

## 6. Proposed Results Outline

- Software architecture and RFU backend integration.
- Single-cell RFU annotation in AnnData.
- Group-level RFU aggregation and visualization.
- Public dataset demonstration.
- Optional radiation PBMC application.

## 7. Target Venues

- APBC short/software paper.
- TCBB extended version.
- Current Bioinformatics.
- Computational Biology and Chemistry.
- Quantitative Biology.

## 8. Immediate Blockers

- Need one real dataset demonstration.
- Need workflow figure.
- Need benchmark or feature table.
- Need internal review-ready manuscript draft.
