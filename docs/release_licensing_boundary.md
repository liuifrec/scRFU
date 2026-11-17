# Software and external-data boundary

The MIT license in this repository applies to scRFU's own code and packaged
synthetic fixtures. It does not replace or reinterpret the terms attached to
external resources.

- The upstream RFU implementation and its reference assets must be obtained
  separately by the user. scRFU records their hashes and invokes them when
  configured, but does not redistribute them.
- VDJdb must be obtained separately from its official release source. Its
  release identifier and database hash belong in runtime provenance; database
  rows are not package fixtures.
- Wells, GEO, and other public validation datasets remain under their source
  repositories' access and usage terms. scRFU distributes neither those inputs
  nor real derived result tables.

Installation tests and the default tutorial use only small synthetic data.
External-resource commands are opt-in and keep outputs outside the repository.
