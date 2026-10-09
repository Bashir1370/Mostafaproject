# Step 05a frozen public metadata reference

These are the exact raw PDBe/SIFTS and AlphaFold DB responses used in the complete implementation validation, acquired on 2026-10-09 UTC. They contain public structure metadata for the accepted 760 canonical mouse proteins, not unpublished experimental data or coordinate files.

reference.json binds the archive and historical source ledger with SHA256. The script validates every pair before import, preserves original retrieval timestamps and runs without any API requests when --reference-cache is selected. Reuse is explicitly marked frozen_reference; it is not a fresh workstation query. Network caches remain separate. See ../../docs/STEP_05A.md for the complete contract and command.

Validation on real frozen inputs: all 48 regression tests passed (including earlier R/QC/DE checks and four reference-mode tests). The complete 760-protein / 1,520-query reference run succeeds with urllib network access forbidden in the test. Corrupt archive hashes, unsafe members and altered existing reference caches are rejected. No actual Ubuntu Step 05a result is claimed by this validation.
