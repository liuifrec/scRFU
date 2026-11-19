# V(D)J API compatibility review

Status: unreleased development proposal for 0.5.0; no release is created.

## Implementation boundary

Five additive functions are exported through `scrfu.pl`; six preparation
functions and six result types are exported through `scrfu.tl`. Drawing lives
in `_vdj_plots.py`, and preparation in `visualization.py`. The supported public
import routes are the existing `pl` and `tl` namespaces. Existing plotting and
analysis function bodies and signatures are unchanged. The 0.4 snapshot is
untouched; the 0.5 snapshot adds only the new surface.

| Drawing API | Preparation result | Quantitative meaning |
| --- | --- | --- |
| `rfu_umap` | `EmbeddingPlotData` | Existing coordinates in source order; display-only selection |
| `rfu_vj_usage` | `VJUsageResult` | Complete joint calls and marginals, including Other/missing mass |
| `rfu_cdr3_spectratype` | `SpectratypeResult` | Actual AA lengths; valid denominators and exclusions |
| `rfu_membership_network` | `MembershipResult` | Confirmed identity-to-RFU edges; complete abundance counts |
| `rfu_longitudinal_presence` | `LongitudinalPresenceResult` | Explicitly ordered donor/sample support; unavailable visits retained |

`ReceptorPlotData` supplies the common cell/receptor/read policy. Preparation
uses exact joins, collision-free identity tuples and explicit selection. Plot
limits cannot alter complete source aggregations. Native AIRR adaptation reuses
the existing canonical adapter and saved-schema checks without assignment or
cell-summary broadcasting. Drawing uses Matplotlib alone and returns Axes.

## Compatibility and scope

- Python 3.10/3.11/3.12 and pandas 2/3 remain in the existing CI matrix.
- Core dependencies, extras and existing scientific defaults are unchanged.
- Matplotlib is imported only when drawing; clean core installations can prepare
  counts and fail with an actionable plotting-extra message if drawing is requested.
- Wheel and source-distribution gates exercise both core-only installation and a
  synthetic gallery after adding the plotting extra, outside the checkout.
- All new public fixtures and gallery content are invented. No participant tables,
  biological embeddings, cohort-specific filters, reference assets or manuscript
  figures are required or included.
- The synthetic tests verify defined counting and rendering behavior, not the
  biological validity of RFU function, antigen recognition, or clone inference.

The ordinary Python suite excludes the existing opt-in external R integration
checks. This plotting change does not require running a reference model. Optional
scverse warnings originate from upstream storage/version interfaces and should
not be interpreted as failed scientific checks.

## Review decisions

The default exact identity is AA CDR3 + literal V call; a caller can explicitly
choose nucleotide/gene tuples or independently called clonotype IDs. Cell mode
rejects ambiguity instead of selecting a primary chain. Read mode requires real
read counts. Missing values are exposed rather than inferred. All category
limits are deterministic and reversible through the returned full tables.

No sequence-distance network backend is included: that requires an independent
distance definition. No mandatory Scirpy, Scanpy, Dandelion, NetworkX or Logomaker
dependency is added. The alluvial is descriptive identity presence, not cell
tracking, causal treatment response or proof of clonal persistence.

An additive **0.5.0** release is appropriate after API review. The checkout uses
**0.5.0.dev0** until that decision; published release tags and assets are immutable.
