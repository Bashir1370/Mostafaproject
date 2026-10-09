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

## Mapping workstation review

Both actual Ubuntu mapping snapshots are uploaded at source commit 3bdf58061337550aeda2f180e44fc710f810170b and reviewed. All snapshot/input hashes, candidate/evidence tables, mapping/FASTA associations and response body hashes pass; combined coverage is 789 genes / 787 proteins. See docs/audits/ubuntu_resolution_2026-10-09/REVIEW.md. The 65 unresolved identities remain open and audited; acceptance is limited to the mapped subset. No coordinate/local-quality assessment, structural score, benchmark or enrichment has run. Original study metadata limitations remain as documented. GSE125002 remains separate and unaudited.

## Step 03 implementation validation

The standard-library Python script retrieved 35 complete UniProt batches from release 2026_03. It audited all 854 genes: 285 pass / 283 unique proteins, 562 held for multiple mouse entries, six held without a current exact link, one unassessable fragment. Cached replay reproduces these counts. All nineteen tests passed, including earlier input/QC/DE checks. See docs/STEP_03.md and docs/audits/Mapping_validation. Held genes are not biological exclusions; the accepted first-pass subset has limited coverage.

## Step 03b implementation validation

Independent gene-centric representative evidence resolves 504 more genes without selecting by review status or length. Combined result: 789 / 854 genes (92.4%) and 787 unique proteins, with 64 held and one fragment unassessable. Baseline mappings are preserved. All twenty-six tests pass, including earlier input/QC/DE and first-pass mapping tests. See docs/STEP_03B.md and docs/audits/Resolution_validation. The relationship service did not supply a release header; response hashes/access times are recorded separately from the UniProt sequence release 2026_03.

## Step 04 implementation validation

Reviewed snapshot inputs are frozen in config/GSE286387_step04_input.json. The offline canonical inventory yields 760 C-positive proteins, 27 no-C proteins, 31 with one C, 729 with multiple C and 10,799 unique sites. The full 854-gene coverage audit retains the 64 held / one unassessable upstream cases. All 31 regression tests passed; the five Step 04 tests also pass in a clean offline fixture without prior results/cache. The actual Ubuntu Step 04 reproduction is now uploaded and reviewed separately at source commit e419a4e98e8b51368a0d5123b08fc128415eee1c; every full sequence/site/audit table and all snapshot/upstream checksums pass. See docs/STEP_04.md and docs/audits/Cysteine_validation.

## Step 04 workstation acceptance

The actual Python 3.13.13 Ubuntu inventory snapshot is reviewed; accepted 760 C-positive proteins / 10,799 sites are bound in config/GSE286387_step05_input.json. Full discovery gene audit retains 760 C-positive mappings, 29 no-C mappings, 64 held and one unassessable. See docs/audits/ubuntu_cysteine_2026-10-09/REVIEW.md. At the user's request, actual results/ outputs are now tracked in Git; reviewed snapshots remain frozen provenance.

## Step 05a implementation

Metadata catalogue queries both PDBe/SIFTS and AlphaFold DB for every accepted protein. All offered chains/models are retained with mouse/accession/sequence gates and raw response SHA256; network errors are not missing structures. Candidate metadata is not structural inclusion: all eligibility fields remain not_assessed and no site passes are issued. Coordinate retrieval, residue mapping, assembly/mutation review and local quality remain Step 05b work. Live 1,520-query validation and final offline replay find 94 proteins with mouse PDBe options, 731 with sequence-matched AlphaFold options and 29 without metadata-eligible candidates. All 42 regression tests pass. See docs/STEP_05A.md and docs/audits/Structure_catalogue_validation.

## Next

The first Ubuntu Step 05a attempt failed with HTTP 403; no final catalogue was generated. Diagnostic revision 1.0.1 records the failing endpoint and bounded response in failure_context.json. Pull and rerun scripts/05a_catalogue_GSE286387.py; send the diagnostic file if it fails, or catalogue_report.md if successful. Root cause of workstation access rejection remains pending. Commit actual successful outputs under results/ after completion. Review the catalogue before Step 05b coordinate and local-quality work. No mapping is forced and no oxidation probability is claimed.
