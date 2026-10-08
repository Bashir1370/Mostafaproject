# Step 04 — canonical cysteine inventory

Protocol 0.1.0; implementation 1.0.0. This is scientific Step 04, not a new scoring stage.

## Input and gates

Use the reviewed combined Ubuntu Step 03b snapshot at source commit 3bdf58061337550aeda2f180e44fc710f810170b, bound by config/GSE286387_step04_input.json. The script verifies all 20 historical snapshot files plus the upstream discovery, configuration and mapping scripts before counting. It runs offline with Python >=3.8 standard library; no R package or network request is needed.

Input coverage: 854 discovery genes, 789 accepted gene mappings, 787 unique canonical mouse accessions. The 64 held / one unassessable mappings remain visible in a separate full-universe gene audit. They are not Cys-inventory input proteins.

Entry requires taxid 10090, displayed canonical sequence, matching sequence length/SHA256, valid uppercase alphabet, one FASTA record per accepted accession, exact FASTA-to-table gene association and the frozen gene universe. Standard amino acids plus U/O are allowed as in Step 03; U is not counted as C. Missing, ambiguous, duplicate or altered inputs stop the run and produce FAILURE.txt rather than a partial success.

## Fixed inventory rule

- Count only residue C; count every occurrence in the unchanged complete displayed sequence.
- Use 1-based canonical numbering; retain initiator, signal-peptide and propeptide positions. No mature-chain renumbering occurs.
- Site identity is accession:Cposition; genes sharing an accession do not duplicate sites.
- Protein passes to prospective Step 05 if it contains at least one C; zero C gives excluded / NO_CYSTEINE_IN_ACCEPTED_SEQUENCE.
- Calculate descriptive Cys percentage as 100*n_cys_total/sequence_length. Neither count nor percentage enters the frozen structural score.
- Do not filter disulfides, metal-binding or modified sites yet; those belong to Step 06. Sequence C does not establish a free thiol, oxidation or structural assessability.

## Ubuntu command

From /home/bashir/Desktop/Mostafaproject:

```bash
git pull --ff-only
python3 oipn-cysteine-structural-screening/scripts/04_inventory_GSE286387.py
cat oipn-cysteine-structural-screening/results/04_cysteine_inventory/inventory_report.md
```

## Outputs

| File | Contract |
|---|---|
| protein_inventory.csv | One row per 787 accepted accessions, including zero-C proteins; length, sequence SHA, all associated stable genes, total C, percentage, positions, pass/excluded and reason. |
| cysteine_inventory.csv | One row per unique C site; site_id, accession, canonical_cys_position, associated genes, total C, percentage, sequence length/SHA, release and protocol. |
| step_audit.csv | One terminal protein status per accepted input accession. |
| gene_step_audit.csv | One row per 854 discovery genes; C-based decision for mapped genes, inherited held/unassessable status for unassessed upstream cases. |
| cysteine_positive_sequences.fasta | Original sequences for C-positive accessions only; gene associations retained. |
| input_checksums.csv | Reviewed binding, all frozen input files and executing Step 04 script SHA256. |
| summary.json; inventory_report.md | Separate gene/protein/site counts and environment provenance. |
| SUCCESS.txt / FAILURE.txt | Success requires STEP04_GENERATED_REVIEW_PENDING and no FAILURE.txt. |

Output is staged and success written last; an unsuccessful rerun removes the old success marker. Older tables may remain after failure and must not be consumed. Historical mapping outputs and configuration are not overwritten.

## Implementation validation, not workstation reproduction

The reviewed workstation sequence input yields 760 C-positive proteins, 27 zero-C proteins, 31 with exactly one C, 729 with multiple C and 10,799 unique sites. This execution occurred in the validation environment; the user's Ubuntu Step 04 output is pending. Tests cover terminal/adjacent positions, U handling, gene-to-accession deduplication, invalid FASTA, changed frozen inputs/failure markers, and all real-sequence positions checked against an independent regex-based enumeration.

Review the Ubuntu report and freeze its accepted protein/site input before Step 05. No structure availability, chemical-state classification, SASA, pKa or structural score has been computed here. No-Cys exclusion applies only to this Cys-based method and says nothing about other oxidation chemistry.
