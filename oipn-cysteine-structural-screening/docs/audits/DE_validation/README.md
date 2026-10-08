# Step 02 implementation validation — 2026-10-08

These are implementation-validation records from the official counts with R 4.3.3 / DESeq2 1.42.0, not the user-workstation DESeq2 1.42.1 result. All eleven tests passed in 37.047 seconds. Independent Python checks recomputed BH over all valid biotypes and verified complete auditing and exact coding significance selection. Modified counts and accepted sample manifest were rejected before fitting. Existing count-integrity/QC integration checks also passed.

Input: 40,481 genes, ten samples; retained: 18,538; valid p-values: 18,503; significant all biotypes: 893; coding: 854 (477 up, 377 down). All 35 missing post-filter p-values were recorded as unassessable rather than silently discarded. Full runtime outputs remain in ignored results. No structural analysis has run. The small summary, input-checksum and environment snapshots document this validation; local reproduction is the next gate.
