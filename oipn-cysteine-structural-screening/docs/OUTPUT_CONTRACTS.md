# Output contracts and provenance

Steps 01–02 now have executable implementations; Steps 03–10 remain output schemas for future work. Implementation-validation outputs and user-workstation outputs are distinguished in STATUS.md. CSV uses UTF-8, one header, explicit missing values and stable IDs. Each output records protocol version; large files retain a checksum/source manifest.

| Step | Main output | Required fields/content |
|---|---|---|
| 01 | sample_manifest.csv; source_manifest.csv; count_matrix.tsv; qc_report.md | sample_id, GEO_sample_id, animal_or_pool_id, condition, species, tissue, dose, route, collection_time, sex, age, strain, batch, biological_unit, inclusion_status, reason; raw-count provenance |
| 02 | all_genes_de.csv; significant_protein_coding_degs.csv; tested_gene_universe.csv | stable_gene_id, symbol, biotype, baseMean, log2FC, lfcSE, statistic, pvalue, padj, prefilter_status, na_reason |
| 03 | gene_protein_mapping.csv; canonical_sequences.fasta | stable_gene_id, protein_accession, isoform_id, canonical_status, species_taxid, sequence_length, sequence_checksum, reference_release, mapping_status, reason |
| 04 | cysteine_inventory.csv | site_id, protein_accession, canonical_cys_position, n_cys_total, cys_percent, sequence_length |
| 05 | structure_manifest.csv; site_structure_mapping.csv | site_id, source, structure_id, checksum, model, assembly, chain, structure_residue_number, insertion_code, SG_present, sequence_identity_context, local_completeness, target_plddt, minimum_neighbor_plddt, quality_status, reason |
| 06 | site_state_audit.csv | site_id, chemical_state, evidence_source, evidence_level, annotation_completeness, primary_eligible, reason |
| 07 | cysteine_features.csv | site_id, SG_SASA_A2, predicted_pKa, nearest_other_SG_distance_A, distance_missing_reason, calculation_status, structure_checksum, tool_versions, preparation_id |
| 08 | site_scores.csv; protein_ranking.csv; shortlist.csv; rank_reference.csv | site_id, R_SASA, R_minus_pKa, site_score; gene/protein IDs, protein_score, tied_best_sites, tied_rank, n_cys_total, n_cys_covered, n_cys_scored, scored_cys_coverage, source, shortlist_threshold, shortlisted |
| 09 | benchmark_audit.csv; performance.csv; robustness.csv; interpretation_status.md | experimental_site_mapping, chemistry, labels/continuous response, measurement_status, metrics/uncertainty, baseline_comparison, scenario, reference_membership, rank_shift, eligibility_change, shortlist_retention |
| 10 | final_candidates.csv; enrichment.csv; final_report.md | rank, gene/protein, best_Cys, structural_score, coverage, confidence, log2FC, padj, location, function, redox_evidence, replication_context, robustness, limitations; enrichment background and tested term collection |

## Mandatory step audit

step_id, object_level (sample/gene/protein/site), object_id, status (pass/excluded/held/unassessable), reason_code, reason_text, input_source, input_checksum, protocol_version. Every input object has one terminal status per step at its object level. Export summary counts separately for sites and proteins. Excluded/unassessable objects remain in the audit.

## Tool/environment provenance

Record Ubuntu release, architecture, R sessionInfo, Python package versions, tool arguments, retrieval date, database release, input/reference checksums and random seeds. A plan to use a tool is not a tested version lock.

## Step 02 details

`all_genes_de.csv` includes every original unique gene, including prefilter exclusions. `tested_gene_universe.csv` is the retained test set before p-value missingness; BH uses its finite p-values. `significant_protein_coding_degs.csv` contains exactly coding padj<0.05 rows. Additional fields are annotation_reference, maxCooks, valid_p_for_BH, significant, passes_to_step03, direction and protocol_version. Missing results have explicit reasons. Supporting outputs include sample_manifest_used.csv, design_matrix.csv, input_checksums.csv, summary.json, de_report.md, sessionInfo.txt and dds_fitted.rds. Runtime success/failure markers govern eligibility for downstream use.
