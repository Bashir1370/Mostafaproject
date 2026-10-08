# Ubuntu Step 02 snapshot — 2026-10-08

Status: STEP02_SNAPSHOT_ACCEPTED_FOR_STEP03.

The real workstation snapshot was uploaded in commit d80a26859ad54577ebf4021e87ef4f4a0cb8cac0. It includes the report, summary, input checksums, session information, 854-gene discovery CSV, sample manifest, model matrix, success marker and snapshot checksum manifest.

All eight file checksums and the selected-gene/model/provenance checks passed. R 4.3.3 / DESeq2 1.42.1 produced 477 up and 377 down coding DEGs; exact discovery membership matches the separate 1.42.0 validation run. See [review evidence and limits](REVIEW.md).

Step 03 uses the uploaded discovery file identified in config/GSE286387_step03_input.json. Full runtime gene tables and fitted model remain on the workstation. No mapping or structural analysis has run.
