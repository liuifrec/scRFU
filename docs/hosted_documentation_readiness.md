# Hosted-documentation readiness

The repository now contains a Read the Docs v2 configuration that installs the
`docs` extra, uses Python 3.11, points to `docs/conf.py`, and fails on warnings.
Local `sphinx-build -E -W --keep-going` succeeds and writes `objects.inv`.

Hosted builds enable intersphinx inventories for AnnData, Scanpy, MuData, and
Scirpy. Offline local builds omit network inventory retrieval while retaining
the extension and producing scRFU's own inventory.

External actions remain:

1. create/link the documentation service project;
2. verify the default and stable-version builds;
3. choose the durable public documentation URL;
4. replace the prospective URL in the ecosystem metadata draft if needed;
5. verify that the published `objects.inv` is reachable.

No external documentation project or deployment was created in this sprint.

