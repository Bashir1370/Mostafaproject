# Progress checkpoint

Protocol: 0.1.0. Audit recorded: 2026-10-08. Workstation target: Ubuntu Linux.

## Completed

- Ten-step design, structural scoring specification, validation plan and Ubuntu instructions.
- Step 01 official GEO metadata and raw count download/integrity audit executed.
- Exact sample-title to GSM mapping: five vehicle, five oxaliplatin samples.
- 40,757 source rows reduced to 40,481 unique stable IDs by retaining one identical count vector per gene; original annotation alternatives and 276 redundant rows are recorded.
- Executable standard-library Python audit, four integrity regression tests, and actual small audit snapshots committed.

## Pending

Step 01 is **DATA_READY_QC_PENDING**, not completed: R sample QC, animal/pooling identity review, age/batch availability and final sample inclusion. No samples have been removed or finally approved. No DE fitting, sequence/structure retrieval, scoring, validation or enrichment has run. R is unavailable in the assistant execution environment; R execution will be on the user's Ubuntu workstation.

## Next

Reproduce the audit locally using README commands. Then prepare/run R QC (library sizes, transformed expression PCA/sample distances and design review). Freeze sample inclusion and model before Step 02 DESeq2. Missing batch is not evidence of absence of batch effects. GSE125002 remains a separately audited complementary candidate.
