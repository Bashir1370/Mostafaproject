# Step 01 data audit

- status: DATA_READY_QC_PENDING
- input_rows: 40757
- unique_genes: 40481
- n_samples: 10
- duplicate_gene_ids: 167
- redundant_rows_removed: 276
- protein_coding_genes: 21534
- all_zero_genes: 11082
- genes_passing_planned_filter: 18538

Duplicate stable IDs are collapsed only when all counts and biotype agree. Counts are never summed or rounded. Original annotation alternatives remain in source_gene_annotations.csv. No sample is excluded or approved before R QC. No DEGs are reported.

- R sample-level QC before inclusion freeze
- animal/pooling identities, age and batch not specified in GEO sample metadata
- reference is Ensembl release 76; version-aware later mapping required
- GEO mentions max-count filtering; export includes genes below 100, so stage of original filtering is unresolved
