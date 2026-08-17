# Deterministic parallel RFU chunks

RFU assignment remains serial by default. Parallel execution requires a
positive `chunk_size`, `max_workers > 1`, and `executor="process"` or
`"thread"`. The caller, not CPU discovery, chooses worker count.

Parallelism occurs only across independent exact-sequence query chunks. Each
chunk has a deterministic identity, isolated input/output/pending/log files, and
validated manifest. Final concatenation is by chunk index rather than completion
order. Worker count and executor are provenance fields but are excluded from the
scientific cache identity, allowing valid completed chunks to be reused across
serial and parallel orchestration.

On failure, the raised error names the chunk; that chunk is marked failed;
completed chunks remain valid; pending outputs are never accepted as complete.
Parallel mode is experimental pending bounded real scaling and is not the
default performance claim.
