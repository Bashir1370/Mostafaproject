# Progress checkpoint

Protocol: 0.1.0. Audit recorded: 2026-10-08. Workstation target: Ubuntu Linux.

## Completed

- Ten-step protocol and structural scoring specification frozen.
- Official GSE286387 inputs audited: 40,481 unique stable IDs; identical redundant annotation rows collapsed without summing.
- User's Ubuntu R 4.3.3 / DESeq2 1.42.1 sample QC reproduced and committed under `docs/audits/ubuntu_2026-10-08`.
- Full study methods and sample QC reviewed: Step 01 accepted with documented metadata limitations. All ten libraries retained; five control and five oxaliplatin. Frozen unpaired model `~ condition`.
- Accepted manifest and checksum-bound design committed in config; see `docs/DATASET_DESIGN_AUDIT.md`.
- Step 02 implemented and run on real inputs in R 4.3.3 / DESeq2 1.42.0: 18,538 tested, 18,503 valid p-values, 854 significant coding genes (477 up / 377 down). This is implementation validation, not the user's final workstation DE run.
- Eleven input-integrity, QC integration, BH/selection-contract and frozen-input gate tests passed.

## Pending

User-workstation Step 02 execution with DESeq2 1.42.1 and output review. No sequence/structure retrieval, structural ranking, validation or enrichment has run. Original animal identifiers, batch, collection age and exact final-dose interval remain limitations, documented in the accepted design. GSE125002 remains a separate unaudited complementary candidate.

## Next

Pull the update, execute `scripts/02_DE_GSE286387.R` in RStudio or Rscript, and review `results/02_differential_expression/de_report.md` and `summary.json`. After workstation reproduction, proceed to Step 03 version-aware mouse gene/protein mapping.
