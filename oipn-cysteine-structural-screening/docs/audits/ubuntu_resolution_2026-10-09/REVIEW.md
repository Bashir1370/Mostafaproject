# Workstation mapping review — accepted subset for Step 04

Source commit: `3bdf58061337550aeda2f180e44fc710f810170b`. Review recorded 2026-10-09 (Asia/Tehran).

Both Ubuntu mapping snapshots are present. All 20 files listed in their two snapshot manifests pass SHA256 verification. Every upstream configuration/discovery/script binding is verified; Step 03b baseline result references are checked against the uploaded first-pass snapshot, not substituted workstation paths.

Both mapping tables, both canonical FASTAs, both terminal audits, the full first-pass candidate table and gene-centric evidence table are byte-identical to independent implementation-validation outputs. All 35 UniProt and 58 relationship response URL/body-SHA pairs match the independently cached, checksum-verified responses. Access times are independently recorded and need not match. UniProt sequence release is 2026_03; the relationship service supplied no release header.

The accepted discovery universe contains 854 distinct stable genes. Baseline has 285 accepted genes / 283 accessions; 504 additional genes pass the documented explicit representative checks. Combined: 789 accepted genes / 787 unique accessions, 64 held and one unassessable. All baseline accepted mappings remain unchanged; every accepted sequence checksum and FASTA gene association is verified. Gene-level RNA does not establish the expressed isoform.

## Remaining cases and disposition

| Reason | Genes | Disposition |
|---|---:|---|
| REPRESENTATIVE_LACKS_ORIGINAL_EXACT_GENE_LINK | 49 | Held; an outside representative needs a separately documented identity bridge. |
| REPRESENTATIVE_ONLY_NONCANONICAL_OR_UNRESOLVED_ISOFORM_XREF | 7 | Held; current displayed-sequence association is not established. |
| NO_CURRENT_EXACT_ENSEMBL_XREF | 6 | Held; historical-ID review remains open. |
| REPRESENTATIVE_MULTIPLE_DISPLAYED_ISOFORMS | 1 | Held; displayed-sequence ambiguity remains. |
| MULTIPLE_GENE_CENTRIC_REPRESENTATIVES | 1 | Held; conflicting representative evidence remains. |
| FRAGMENT_SEQUENCE | 1 | Unassessable; incomplete sequence is not repaired. |

This is review of the evidence gates and acceptance of the 789-gene subset, not manual rescue or closure of the other 65 identities. No missing case is counted as Cys-free or oxidation-resistant. The requested majority coverage is achieved descriptively (92.4%); no mapping is forced to reach a threshold.

The accepted mapping/FASTA and snapshot/provenance hashes are frozen in `config/GSE286387_step04_input.json`. Step 04 uses the committed reviewed snapshot offline; raw network caches are not required. Generation markers in the historical snapshots remain unchanged. Structure retrieval and susceptibility ranking are not approved by this sequence-level review.
