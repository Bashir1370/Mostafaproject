# Step 03b implementation validation — 2026-10-08

The accepted 854-gene discovery input, original Step 03 outputs and checksum-verified UniProt release 2026_03 entries were used. Fifty-eight complete accession query batches were retrieved from the EMBL-EBI gene-centric service. Its release header was absent; record response hashes/access times separately, without claiming a matching relationship release.

504 additional genes were resolved by explicit gene-specific grouping, one representative, complete candidate accounting and original mouse/canonical/sequence checks. Combined result: 789 genes (92.4%), 787 proteins, 64 held and one unassessable fragment. No reviewed/longest-entry shortcut or symbol rescue was applied. Original 285 accepted gene mappings are unchanged. Full cached replay reproduced totals.

All 26 tests passed in 38.328 seconds, including earlier Steps 01–03 tests and seven resolution/baseline-integrity/output tests. The final resolver/helper/input hashes are recorded in input_checksums.csv. These are implementation-validation artifacts, not the user's Step 03b workstation files. Workstation reproduction and review are pending; no Cys filtering or structural scoring has run.
