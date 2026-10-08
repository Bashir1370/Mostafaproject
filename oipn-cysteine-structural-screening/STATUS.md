# Progress checkpoint

Protocol: 0.1.0. Audit recorded: 2026-10-09 (Asia/Tehran). Workstation target: Ubuntu Linux.

## Completed

- Ten-step protocol and structural scoring specification frozen.
- Official GSE286387 inputs audited: 40,481 unique stable IDs; identical redundant annotation rows collapsed without summing.
- User's Ubuntu R 4.3.3 / DESeq2 1.42.1 sample QC reproduced and committed under `docs/audits/ubuntu_2026-10-08`.
- Full study methods and sample QC reviewed: Step 01 accepted with documented metadata limitations. All ten libraries retained; five control and five oxaliplatin. Frozen unpaired model `~ condition`.
- Accepted manifest and checksum-bound design committed in config; see `docs/DATASET_DESIGN_AUDIT.md`.
- Step 02 implemented and run on real inputs in R 4.3.3 / DESeq2 1.42.0: 18,538 tested, 18,503 valid p-values, 854 significant coding genes (477 up / 377 down). This is implementation validation, not the user's final workstation DE run.
- Eleven input-integrity, QC integration, BH/selection-contract and frozen-input gate tests passed.

## Workstation reproduction and review completed

The actual Ubuntu R 4.3.3 / DESeq2 1.42.1 snapshot is uploaded and reviewed at source commit d80a26859ad54577ebf4021e87ef4f4a0cb8cac0. All eight snapshot checksums, configuration provenance, sample identities and model matrix pass. Discovery has 854 unique significant coding genes (477 up / 377 down); exact membership matches implementation validation. See docs/audits/ubuntu_DE_2026-10-08/REVIEW.md. The accepted discovery file/checksum is bound in config/GSE286387_step03_input.json.

## Pending

Both actual Ubuntu mapping snapshots are uploaded at source commit 3bdf58061337550aeda2f180e44fc710f810170b and reviewed. All snapshot/input hashes, candidate/evidence tables, mapping/FASTA associations and response body hashes pass; combined coverage is 789 genes / 787 proteins. See docs/audits/ubuntu_resolution_2026-10-09/REVIEW.md. The 65 unresolved identities remain open and audited; acceptance is limited to the mapped subset. No structure, structural score, benchmark or enrichment has run. Original study metadata limitations remain as documented. GSE125002 remains separate and unaudited.

## Step 03 implementation validation

The standard-library Python script retrieved 35 complete UniProt batches from release 2026_03. It audited all 854 genes: 285 pass / 283 unique proteins, 562 held for multiple mouse entries, six held without a current exact link, one unassessable fragment. Cached replay reproduces these counts. All nineteen tests passed, including earlier input/QC/DE checks. See docs/STEP_03.md and docs/audits/Mapping_validation. Held genes are not biological exclusions; the accepted first-pass subset has limited coverage.

## Step 03b implementation validation

Independent gene-centric representative evidence resolves 504 more genes without selecting by review status or length. Combined result: 789 / 854 genes (92.4%) and 787 unique proteins, with 64 held and one fragment unassessable. Baseline mappings are preserved. All twenty-six tests pass, including earlier input/QC/DE and first-pass mapping tests. See docs/STEP_03B.md and docs/audits/Resolution_validation. The relationship service did not supply a release header; response hashes/access times are recorded separately from the UniProt sequence release 2026_03.

## Step 04 implementation validation

Reviewed snapshot inputs are frozen in config/GSE286387_step04_input.json. The offline canonical inventory yields 760 C-positive proteins, 27 no-C proteins, 31 with one C, 729 with multiple C and 10,799 unique sites. The full 854-gene coverage audit retains the 64 held / one unassessable upstream cases. All 31 regression tests passed; the five Step 04 tests also pass in a clean offline fixture without prior results/cache. Implementation validation is distinct from the pending Ubuntu Step 04 reproduction. See docs/STEP_04.md and docs/audits/Cysteine_validation.

## Next

Pull and execute scripts/04_inventory_GSE286387.py on Ubuntu; send inventory_report.md for review. Commit/review the workstation inventory before freezing Step 05 input. Structure-retrieval implementation is still pending. No mapping identities are forced and no oxidation probability is claimed.
