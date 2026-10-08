# Step 02 workstation snapshot review — 2026-10-08

Status: STEP02_SNAPSHOT_ACCEPTED_FOR_STEP03
Reviewed source commit: d80a26859ad54577ebf4021e87ef4f4a0cb8cac0.

## Verified

- All eight snapshot file SHA256 values match snapshot_checksums.sha256; SUCCESS.txt records STEP02_COMPLETE.
- Configuration, frozen design and accepted sample-manifest hashes match the workstation input-checksum record. The reported count-matrix checksum matches the frozen design and implementation-validation input.
- Sample records reproduce all ten accepted samples, identities and conditions. The model matrix is an intercept and oxaliplatin indicator; control remains the reference, five samples per group.
- sessionInfo.txt confirms Ubuntu 24.04.4 LTS, R 4.3.3 and DESeq2 1.42.1.
- Discovery CSV contains 854 unique stable gene IDs, each historical protein_coding, prefilter-passing and finite BH padj <0.05. There are 477 up and 377 down; signs agree with the contrast. Rows are sorted by padj and stable ID.
- All 854 selected stable IDs match the separately validated real-data DESeq2 1.42.0 run. Numeric values differ slightly between environments; matching membership does not imply byte-identical numeric results.
- 700 selected genes have absolute log2FC below 1, consistent with the registered absence of an absolute fold-change cutoff.

## Accepted downstream input

Use significant_protein_coding_degs.csv in this snapshot as the Step 03 discovery input, SHA256 c9bf672d38f877361b3132065884d3a02e3132bcbf9901656b9f2cb04fdf3130. The input binding is recorded in config/GSE286387_step03_input.json. Gene IDs, rather than symbols, define objects.

## Review limits and next step

The full all-gene output, tested universe, per-gene audit and fitted RDS remain on the workstation; they were not uploaded in this selected snapshot. This review did not recompute workstation BH from the complete p-value universe. BH was independently verified in implementation validation, and the uploaded discovery membership agrees with that run. The study's public-metadata limitations remain unchanged.

Step 03 will map these genes to mouse canonical proteins and sequences with explicit reference/version provenance. No protein count, mapping result or sequence retrieval is claimed by this acceptance. Unmapped or ambiguous genes must be recorded rather than silently dropped; 854 genes need not become 854 proteins.
