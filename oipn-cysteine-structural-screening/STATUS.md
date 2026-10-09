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

## Step 05a workstation acceptance

Actual Ubuntu catalogue uploaded at f8a1a4e8c60593101425ebc09ba55c2081393bba is reviewed. All eight output files, 42 input checksums and four complete data/audit tables pass; tables reproduce the frozen public reference exactly. Accepted snapshot is docs/audits/ubuntu_catalogue_2026-10-09 and prospective Step 05b binding is config/GSE286387_step05b_input.json. Scope: 760 proteins / 10,799 prospective sites; 731 proteins with sequence-matched AlphaFold options, including 94 with mouse PDBe options; 29 without metadata-eligible options. 27 are unassessable; two retain unresolved offered metadata and are held. Metadata reuse is explicitly frozen_reference, not fresh workstation acquisition. No structure/site is approved.

## Step 05b1 acquisition implementation

Resumable downloader verifies the entire accepted snapshot/upstream input chain and plans 1,176 unique coordinate files (445 PDB entries / 731 AlphaFold models) plus 1,462 offered AlphaFold confidence/PAE files. Every metadata record remains in the audit. Real four-kind probe and acquisition/cache/error/archive tests pass; the bulk workstation download is pending. Archive read-back checks every stored member. Raw acquisition does not approve assemblies, sequences, chemical state or local Cys quality. See docs/STEP_05B1.md and docs/audits/Structure_download_validation.

## Step 05b1 actual incomplete acquisition

Uploaded source ed5ae356e169d2ae5be5b657bbce5130de582bc6 reconciles 2,638 unique requests: 2,121 successful, 517 failed. All 731 AlphaFold coordinates and all 731 confidence files succeeded; 438 PDB coordinates / 221 PAE files succeeded. Remaining: six PDB SSL handshake timeouts, one interrupted large PDB read, 498 PAE HTTP 400 and 12 PAE HTTP 502. Original HTTP 400 bodies were not captured, so origin/root cause remains unresolved. The uploaded receipt references a partial local archive, not an independently inspected remote tar.gz. See docs/audits/ubuntu_download_incomplete_2026-10-09/REVIEW.md.

## Step 05b1 actual complete acquisition accepted

Source b2205ac65ebf0881653dd985ee279aef2f9fa994 is reviewed: all 2,638 files have successful acquisition status, including all 1,176 coordinates. All 53 input checksums and complete plan/candidate associations pass. The previous 2,121 successful raw files retain identical acquisition hashes/metadata. Every raw/sidecar hash reconciles with the 5,337-row archive manifest and 5,338-member receipt. Actual local tar.gz/body bytes were not remotely inspected. See docs/audits/ubuntu_download_2026-10-09/REVIEW.md; accepted immutable manifests and a mandatory local hash gate are bound in config/GSE286387_step05b2_input.json. Earlier partial-run evidence remains historical provenance, not current incompleteness.

## Step 05b2 mapping implementation

Offline Gemmi 0.7.5 mapper verifies every accepted raw file/sidecar before parsing, separates canonical/label/author numbering, checks exact predicted sequences or strict deposited/unique-exact reference correspondence, and preserves every candidate model/subchain. Expected ambiguity and unavailable SG/residues are explicitly audited. All 760 proteins / 10,799 sites remain represented; metadata held/excluded records remain in the 1,872-row candidate audit. Validation: 73 full regression tests and all 14 mapping tests (including an additional end-to-end representative integration) pass. Real accepted-byte examples include 11gl chains A/B and AF-A0A087WRH0-F1; the complete 1,176-file workstation run remains pending. See docs/STEP_05B2.md and docs/audits/Structure_mapping_validation.

## Step 05b2 workstation mapping accepted

Actual source 7dd9d5dde4e2a887a3af46f5c6e0c7b90889641c is reviewed. All nine output Git hashes and complete frozen input checksums pass. Reconstructed 10,799-site audit matches exactly; every candidate/raw association and 2,638 raw/sidecar ledger hashes reconcile. 9,741 unique sites in 731 proteins have mapped SG options (all have AlphaFold; 394 additionally have experimental options). 265 sites remain mapping-unresolved; 793 lack metadata-eligible structures. See docs/audits/ubuntu_cysteine_mapping_2026-10-09/REVIEW.md. Raw byte verification/full parsing occurred on Ubuntu; full raw bodies remain local. This accepts mapping only.

## Step 05b3 local-confidence implementation

An offline local-confidence/target-completeness gate now uses the frozen >=90 target / >=70 modeled-neighbor within 6 A criteria, checks JSON/CIF confidence, retains every mapped experimental option as held for manual local validation, and audits all canonical sites. PAE is descriptive with no new cutoff; neighbor >=90 remains a sensitivity flag. Local pass does not approve assemblies, native/fragment context, neighbor completeness or chemical state. Validation: 90 full regression tests plus a final rerun of all 16 local-quality tests; real 35-site sample and full-universe representative integration pass. See docs/STEP_05B3.md.

## Next

Run scripts/05b3_local_quality_GSE286387.py with the existing subproject .venv Python and upload results/05b3_local_quality. Review actual model/site confidence distributions and holds, then continue experimental local validation, assembly/native/mutation/fragment and neighbor-completeness context. No final Step 05 pass, features or score has been issued.
