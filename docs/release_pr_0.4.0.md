# Proposed release PR

## Title

Release scRFU 0.4.0

## Body

### Scope

This PR prepares scRFU 0.4.0 as the first coherent release after the historical
v0.2.0 and v0.3.0 tags. It promotes package metadata from 0.4.0rc1 to 0.4.0 and
finalizes the changelog, release notes, installation guidance, citation
metadata, artifact boundary, and release scorecard. It does not change frozen
TCR RFU scientific definitions.

### Validation

- 349 tests collected: 345 passed and 4 optional integrations skipped in the
  reviewed release-candidate baseline.
- Exact official-RFU parity and Linux integration are validated; external RFU
  assets remain user supplied.
- GitHub Actions run 32807983086 passed Linux Python 3.10/3.11/3.12, Ubuntu,
  macOS, and Windows Python-only smoke, documentation, optional dependencies,
  and package installation for the compatibility commit.
- Wheel and sdist clean-install, synthetic tutorial, warnings-as-errors Sphinx,
  artifact-content, and reproducible-build checks pass on the prepared local
  tree and must remain green on the final release commit.

### External dependency boundary

The package redistributes no upstream RFU implementation/assets, RFU `.Rdata`,
VDJdb database, Wells/GEO/OAS/10x data, or real-data outputs. Official RFU/R
execution is validated on the primary Linux environment. Cross-platform CI
covers Python-only package functionality and does not claim RFU/R execution on
macOS or Windows.

### BCR scope

Experimental BCR canonicalization, QC, pairing, and feature extraction remain
available. The BCR functional-unit feasibility gate is NO-GO; no BCR reference
or assignment API is included.

### Known limitations

- Real RFU assignment requires an external compatible RFU checkout and R.
- VDJdb analyses require a separately obtained, version-pinned local reference.
- Threshold-qualified assignment is a frozen-reference similarity policy, not
  calibrated confidence.

### Release checklist

- [ ] Maintainer reviews the complete release diff.
- [x] Runtime, wheel, sdist, CITATION, changelog, and release notes report 0.4.0.
- [x] Full tests, Ruff, diff check, build, and warnings-as-errors docs pass.
- [x] Official and development RFU Linux integrations pass with expected skip.
- [x] Clean wheel and sdist tutorial/smoke checks pass outside the source tree.
- [x] SOURCE_DATE_EPOCH builds are byte-for-byte reproducible.
- [x] Wheel and sdist content/privacy audits pass.
- [ ] Remote CI passes on the exact release commit.
- [ ] No merge, tag, upload, release, archive, or DOI occurs without separate
      authorization.
