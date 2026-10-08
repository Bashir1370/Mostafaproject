# Ubuntu Step 02 snapshot — 2026-10-08

## Current evidence

The user supplied the text of de_report.md after running the registered Step 02 script on Ubuntu with R 4.3.3 / DESeq2 1.42.1. This report is recorded verbatim. Its totals match the separate implementation-validation run: 18,538 tested, 18,503 valid p-values, 893 significant genes, 854 significant coding genes (477 up / 377 down).

## Workstation file transfer and review

At creation, only the user-supplied report and this receipt are available in this directory. Upload the actual summary.json, input_checksums.csv, sessionInfo.txt, significant_protein_coding_degs.csv, sample_manifest_used.csv, design_matrix.csv and SUCCESS.txt from the workstation using [the documented commands](../../STEP_02.md). The commands also copy the actual report and generate snapshot_checksums.sha256. Their appearance in the directory confirms transfer; acceptance still requires review of contents and provenance.

No implementation-validation file has been substituted for a workstation file. Matching totals alone do not verify identical selected genes. Full all-gene tables, gene audit, tested universe and fitted model remain in the workstation's ignored results directory. The original protocol and scripts reproduce them. Step 03 is pending discovery-list/provenance review.
