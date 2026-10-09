# Output contracts and provenance

Steps 01–04 and the Step 05a metadata catalogue have executable implementations; raw coordinate acquisition Step 05b1 is implemented; canonical residue mapping Step 05b2 is implemented; local-quality Step 05b and Steps 06–10 remain future implementations. Implementation-validation outputs and user-workstation outputs are distinguished in STATUS.md. CSV uses UTF-8, one header, explicit missing values and stable IDs. Each output records protocol version; large files retain a checksum/source manifest.

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

## Step 03 details

One gene_protein_mapping.csv row and one step_audit.csv row per discovery gene. Pass/held/unassessable are explicit; no artificial zero susceptibility. mapping_candidates.csv retains every exact mouse-entry option including ambiguous alternatives. Accepted columns include canonical_gene_association, entry_type, entry_version, sequence_version, versioned Ensembl gene/transcript/protein IDs and source/deposited reference distinctions. canonical_sequences.fasta deduplicates accepted accessions while retaining all associated gene IDs. sequence_checksum is SHA256 of the uppercase unwrapped amino-acid string, not FASTA bytes. source_manifest.csv records every completed query page, response checksum, UTC access time and release/date. input_checksums.csv includes the executable script. Runtime SUCCESS.txt is STEP03_GENERATED_REVIEW_PENDING, not final inclusion approval.

## Step 03b combined output

results/03b_mapping_resolution/gene_protein_mapping.csv retains all Step 03 fields plus step03_original_status, resolution_method and resolution_version. Every discovery gene has one terminal step_audit row; all accepted accessions occur exactly once in combined canonical_sequences.fasta. gene_centric_evidence.csv records stable_gene_id, representative_identifier, explicit gene-specific member_accessions, covers_all_original_candidates and resolution_reason. Baseline accepted identities/sequences stay unchanged. source_manifest.csv records relationship response provenance separately from the sequence reference release; absent service release headers are explicit. Success remains generation/review-pending, not final inclusion approval. The combined output becomes the prospective Cys-inventory input only after reproduction/evidence review.

## Step 04 details

The reviewed mapping snapshot is checksum-bound in config/GSE286387_step04_input.json. protein_inventory.csv and protein step_audit.csv include every accepted input accession, including zero-C exclusions. cysteine_inventory.csv deduplicates accession/1-based canonical position, retains gene associations, sequence hash/release and descriptive count/percentage. gene_step_audit.csv retains all discovery genes and inherited upstream held/unassessable cases. Cys-positive FASTA preserves the full displayed sequences. Provenance, summary/report and generation markers follow STEP_04.md. No-Cys exclusion is method-specific; count/density are not score terms.

## Step 05a provisional catalogue

structure_candidates.csv is a service-metadata inventory, distinct from the final structure_manifest.csv/site_structure_mapping.csv after coordinate validation. Its structural_eligibility is always not_assessed. protein_structure_inventory.csv and step_audit.csv cover all 760 accepted accessions, while source_manifest.csv covers both queries per accession. Pending candidates and unresolved offered records are held; no eligible record without an unresolved offered option is unassessable, never zero susceptibility. No site-level structural decision is emitted. See STEP_05A.md for current/legacy field handling and provenance/error gates.

## Step 05b input binding

config/GSE286387_step05b_input.json records STEP05A_SNAPSHOT_ACCEPTED_FOR_STEP05B and binds all actual Ubuntu catalogue snapshot files, their checksum manifest and the Step 05 input binding. Future coordinate scripts must verify that chain and use the frozen candidate/source/protein files, while retaining the canonical 760-protein / 10,799-site universe from the prior binding. The catalogue has no structural eligibility approvals. Step 05b1 acquisition is implemented; subsequent residue mapping/local quality remains pending.

## Step 05b1 acquisition outputs

download_plan.csv has one row per kind/URL, SHA256-derived file_id and required-coordinate flag. candidate_file_links.csv preserves every offered row and its file associations or held/excluded metadata status, always structural_eligibility=not_assessed. download_manifest.csv adds acquisition status/reason, HTTP code, bytes/raw SHA256, original retrieval time, response URL, ETag/Last-Modified and project-relative body/sidecar paths. Network failure is failed; 404 is unavailable, never oxidative resistance. Raw archive_manifest.csv and archive_receipt.json record member/raw/archive hashes separately; partial archives retain incomplete status. PROBE_SUCCESS.txt is not full-download SUCCESS.txt.

Revision 1.0.1 download_manifest.csv additionally retains error_headers (allowlisted header JSON string) and error_excerpt (at most 2,048 response bytes) for HTTP failures, together with explicit HTTP status and response URL. Cookies are omitted; HTTP 400 remains failed. Timeout/attempt and retry_probe controls are recorded in summary.json; the frozen URL/file IDs and raw-cache namespace are unchanged.

## Step 05b2 accepted raw-file input

config/GSE286387_step05b2_input.json binds the complete actual Ubuntu download snapshot, archive receipt/manifest, frozen file-plan SHA256, local raw-file directory and preceding Step 05b binding. Status STEP05B1_MANIFEST_ACCEPTED_FOR_STEP05B2 accepts manifest/provenance consistency only. The mapper verifies each raw body and metadata sidecar against the accepted acquisition/archive hashes before parsing. It must preserve the entire canonical protein/site universe and distinguish absent/ambiguous mapped residues from acquisition/identity errors. Step 05b2 mapping is implemented; no local-quality/site eligibility approval exists yet.

## Step 05b2 mapping outputs

site_structure_mapping.csv retains canonical site ID/position, candidate row/source, author and label chains, deposited model, label/auth residue IDs/insertion code, mapping method/status/reason, deposited/observed components, SG existence/alternate/occupancy/coordinates/raw B factor, atom-name gaps, declared sequence-change positions and frozen raw-file ID/hash/path. All structural_eligibility and local_quality_status remain not_assessed. candidate_mapping_audit.csv includes every offered row; site_step_audit.csv and protein_mapping_coverage.csv retain the full canonical universe, including metadata ambiguities and no-candidate cases. raw_file_checksums.csv records all verified raw body/sidecar hashes. Input/environment/summary/report/marker contracts distinguish correspondence generation from downstream structural acceptance.

## Step 05b3 local confidence

site_local_quality.csv retains every mapped SG option, with local_gate_status (pass/excluded/held), reason, target/neighbor pLDDT, neighbor identities/counts, coherent target atom names, confidence/PAE source hashes, descriptive bidirectional PAE and fragment/context flags. structural_eligibility remains held. Experimental options are held for local validation/context; B factors are not converted to pLDDT. The separate site/protein audits retain all original canonical objects and inherited mapping uncertainty. Input/provenance, raw verification and generation markers follow STEP_05B3.md. No final site inclusion or feature selection is emitted.

Implementation 1.0.1 site_step_audit.csv also retains unresolved_mapping_options and unresolved_metadata_review from the full original candidate/mapping universe, not only mapped options. An excluded predicted option cannot terminally exclude a site with unresolved offered evidence. summary.json and quality_report.md report terminal canonical-site counts separately from per-option local-gate counts.

## Prospective Step 05b4 context-review input

config/GSE286387_step05b4_input.json binds all nine actual corrected Ubuntu Step 05b3 files, checksum manifest/review and the prior input binding. Status STEP05B3_LOCAL_CONFIDENCE_ACCEPTED_FOR_CONTEXT_REVIEW is limited to local-confidence generation/provenance. Future scripts must verify this chain, preserve all 760 proteins / 10,799 sites, and retain unresolved offered evidence. No feature-ready structural pass is bound.
