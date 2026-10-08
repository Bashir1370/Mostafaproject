# Step 03 implementation validation — 2026-10-08

The complete registered Step 03 script was run on the accepted 854-gene workstation discovery input against the live UniProt API. All 35 query batches completed, pagination totals agreed, and all pages identified release 2026_03 (02-September-2026). Cached response replay reproduced the same mapping counts. This snapshot documents implementation validation, not workstation Step 03 execution.

Passed: 285 genes / 283 unique proteins. Held: 562 multi-entry genes and six without current exact links. Unassessable: one fragment. All 854 genes have terminal audit rows; FASTA accession membership, sequence lengths and SHA256 values match accepted mapping rows. No Cys presence filter was applied.

The full nineteen-test suite passed (36.032 seconds), including existing Steps 01–02 tests and eight mapping identity/provenance/output checks. The final executable's SHA256 is included in input_checksums.csv. Full raw API responses and full mapping outputs remain in ignored local data/results. Methods and the important held-cohort coverage limitation are in docs/STEP_03.md. No structures or scores were generated.
