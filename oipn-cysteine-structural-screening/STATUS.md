# Progress checkpoint

Protocol: 0.1.0. Audit recorded: 2026-10-08. Workstation target: Ubuntu Linux.

## Completed

- Ten-step design, structural scoring specification, validation plan and Ubuntu instructions.
- Step 01 official GEO metadata and raw count download/integrity audit executed.
- Exact sample-title to GSM mapping: five vehicle, five oxaliplatin samples.
- 40,757 source rows reduced to 40,481 unique stable IDs by retaining one identical count vector per gene; original annotation alternatives and 276 redundant rows are recorded.
- Executable standard-library Python audit, four integrity regression tests, and actual small audit snapshots committed.

## Pending

Step 01 is **DATA_READY_QC_PENDING**, not completed: R sample QC, animal/pooling identity review, age/batch availability and final sample inclusion. No samples have been removed or finally approved. No DE fitting, sequence/structure retrieval, scoring, validation or enrichment has run. User reports R 4.3.3, Bioconductor 3.18 and DESeq2 1.42.1 installed successfully. The QC script has now been exercised independently on actual official counts with R 4.3.3 and DESeq2 1.42.0; user-workstation reproduction is pending.

## Next

Reproduce the audit locally using README commands. Run the committed R QC script (library sizes, transformed expression PCA/sample distances and design review). Freeze sample inclusion and model before Step 02 DESeq2. Missing batch is not evidence of absence of batch effects. GSE125002 remains a separately audited complementary candidate.

## Step 01b implementation

Blind VST sample QC now implemented and run successfully on official inputs: 10 samples and 18,538 retained genes. PCA, library totals, sample distances/correlations and reproducibility records are generated. Seven input-integrity/integration tests pass. No hypothesis tests or sample inclusion freeze were performed. The R QC substep is part of Step 01, not an eleventh scientific step.
