# ALFU Design Note

## Purpose

This note documents a possible future B-cell extension concept. It is design
only and is not implemented in scRFU.

ALFU means Antibody-Lineage Functional Unit.

## Motivation

BCR repertoires contain lineage evolution information through somatic
hypermutation, class switching, heavy/light pairing, and transcriptional
B-cell states. A future ALFU framework could connect antibody sequence lineage
features with AnnData-based single-cell immune profiling.

## Why RFU Cannot Simply Be Reused

RFU is TCR-specific and trained or defined around TCR sequence embeddings and an
RFU atlas. BCR biology has different sequence constraints, lineage structures,
and metadata requirements. Reusing the RFU atlas directly for BCR analysis would
not be scientifically appropriate without a separate method and validation.

## Possible Future ALFU Components

- Heavy-chain CDR3.
- Light-chain pairing.
- V/J usage.
- SHM burden.
- Isotype.
- Clonotype or lineage tree.
- B-cell transcriptional state.
- Memory, plasma, or naive state.

## Proposed Architecture

A future implementation could live under `scrfu.bcr` or `scrfu.alfu` as an
experimental module. It should be clearly separated from the current TCR RFU
workflow and should expose explicit provenance for any external model, atlas, or
lineage-calling backend.

## Current Status

ALFU is design only. It is not implemented and should not be claimed in the
manuscript except as a clearly marked future direction.
