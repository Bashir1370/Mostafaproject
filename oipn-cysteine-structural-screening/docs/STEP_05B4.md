# Step 05b4 — structural context and official validation evidence

Protocol 0.1.0 and all scientific cutoffs/scoring weights are unchanged. This substep collects evidence needed for context and experimental local validation. It issues no structural pass or feature-ready selection. Assemblies are inventoried as deposited recipes, not generated/selected yet.

## Input and Ubuntu execution

config/GSE286387_step05b4_input.json binds the corrected actual Ubuntu 05b3 snapshot (source 1c850c52db0f8e9523e65073ce30496fb48635ce), all its hashes and the previous input chain. Verify every accepted input/executable and reconstruct the full corrected site audit before requests. Recheck all 2,638 raw bodies/sidecars against accepted hashes before coordinate parsing. Keep all original caches. Use existing Gemmi 0.7.5; no new package is required.

```bash
cd /home/bashir/Desktop/Mostafaproject
git pull --ff-only

oipn-cysteine-structural-screening/.venv/bin/python oipn-cysteine-structural-screening/scripts/05b4_context_GSE286387.py --probe-validation &&
oipn-cysteine-structural-screening/.venv/bin/python oipn-cysteine-structural-screening/scripts/05b4_context_GSE286387.py --download-validation --workers 2 &&
cat oipn-cysteine-structural-screening/results/05b4_structure_context/context_report.md
```

The probe tests/caches one official XML report (11gl) in a separate results/05b4_structure_context/probe directory. PROBE_SUCCESS.txt is not full context success. The full execution uses its verified cache and requests the remaining reports. Default timeout 120 seconds / four attempts; --timeout (1..600), --attempts (1..8), --workers (1..4) are transport controls, not scientific cutoffs. Successful cached reports retain original body hashes, HTTP/source URLs, ETag/Last-Modified and retrieval times; no refresh or TLS verification bypass occurs. Interrupted/transient 429/5xx/network reads retry with bounded backoff. 404 is unavailable report evidence; other HTTP/network or identity/hash/parser errors fail execution. Failure diagnostics are in validation_download_manifest.csv and failure_context.json. Rerunning reuses only verified successful pairs.

Running without --download-validation is cache-only and explicitly records uncached reports as not_downloaded; context inventory can be generated but this does not complete experimental validation. Missing/404 reports never indicate absent protein or oxidative resistance.

## Scope and planned sources

Retain all 760 proteins / 10,799 sites, all 14,310 mapped options and the inherited corrected status counts. Parse every one of the 1,176 accepted coordinate files to inventory deposited context. The 445 experimental entry IDs plan 445 official wwPDB XML.gz reports, one per unique entry (not per chain/model). Exact URLs are https://files.wwpdb.org/pub/pdb/validation_reports/<id[1:3]>/<id>/<id>_validation.xml.gz. Compressed raw reports/metadata are stored in data/raw/experimental_validation/step05b4/<plan SHA256>/, outside Git. A successful download requires gzip/XML parsing and exact Entry pdbid; network failure cannot silently become missing structure.

Of the mapped options, inspect the 4,580 predicted local-pass and 4,569 experimental pending options. Retain all 5,161 predicted local-fail rows as unassessed for context, without reintroducing them as eligible. All offered candidate/mapping evidence and its unresolved status remain bound upstream. No held metadata is rescued by a symbol/homolog shortcut.

## Evidence collected, without new eligibility cutoffs

- Preserve entity/source annotations, every assembly definition/generator/operator, declared database sequence differences, deposited connections and unobserved atom/residue annotations. Preserve operator expressions exactly, including products/ranges; do not choose the first assembly or assume its physical expansion was generated. Candidate assembly IDs use label-subchain membership in explicit generator lists.
- For mapped SG, verify exact frozen atom/model/subchain/component/coordinates. Coherent target CYS N/CA/C/O/CB/SG names use only blank atoms plus the selected SG alternate; a blank SG cannot silently select a nonblank side-chain conformer. Do not union different alternate conformers or repair atoms.
- Record all observed non-H/D, nonzero-occupancy atoms within inclusive 6 A of SG in its deposited model, excluding the target residue. Preserve all partner subchains, ligands and alternate states for review. Experimental neighborhoods are deposited-ASU evidence, not biological-assembly neighborhoods. Unknown coordinates prevent claiming a complete geometric neighborhood.
- Assess heavy-atom name gaps separately for each observed alternate of the 20 standard amino acids and SEC; export unrecognized components for manual review. These are descriptive name checks, not geometric or energetic validation. Entirely unobserved residues cannot be located from this calculation; deposited unobserved annotations remain explicit evidence. No missing atom is synthesized.
- Map explicitly documented sequence differences back to canonical positions before flagging local mutations. Author numbering is never substituted for canonical/label numbering. Preserve raw entity mutation annotations and partner species/source annotations; their absence does not prove native state.
- Carry inherited predicted PAE and canonical/fragment flags forward. No PAE threshold, global resolution cutoff or new experimental RSRZ/RSCC/rotamer cutoff is introduced.
- Official ModelledSubgroup records match exact deposited model, author chain/residue, insertion code, component and compatible alternate. Multiple compatible records are ambiguous; absent records/metrics remain missing, not zero. Export raw target/neighborhood attributes and all reported child flags for manual local validation. Do not substitute nearest residue/model or human/rat structure.
- Compare report and coordinate revision dates. Different/unknown dates require compatibility review; equal dates alone do not prove the report validates identical raw bytes. Preserve full Entry/schema attributes and raw report provenance; no global report release is invented.

## Outputs and meaning of completion

results/05b4_structure_context contains validation_download_plan.csv, validation_download_manifest.csv, failure_context.json, site_context_evidence.csv, structure_context_inventory.csv, context_annotations.csv, experimental_validation_residues.csv, validation_entry_inventory.csv, manual_review_queue.csv, site_step_audit.csv, raw_file_checksums.csv, input_checksums.csv, summary.json, context_report.md and generation markers. All dataset site statuses/uncertainty remain inherited, with held sites requiring context review; no final pass is issued. The manual queue groups candidate/subchain/model and counts sites requiring review; pending rows are not approvals and are not automatically consumed for feature extraction. Original raw files are unmodified and preserved. XML reports and their sidecars remain reusable local provenance.

SUCCESS.txt means STEP05B4_GENERATED_REVIEW_PENDING: evidence inventory generated, not all contexts approved. FAILURE.txt invalidates any old tables/SUCCESS. Review actual outputs before binding a subsequent assembly/native/experimental eligibility decision. Suitable experimental mouse retains priority after it passes local/context gates; no structure is selected by final score. Chemical-state classification remains Step 06. No SASA, pKa or oxidation score is generated.

## Validation and practical limits

109 full regression tests passed in 51.600 seconds; the 17 context/report tests passed again in 9.889 seconds after final refinements. Cases cover conformer name gaps, partner/ligand preservation, unknown coordinates, assembly-expression retention, canonical mutation offsets, frozen input/source gates, wrong report identities, corrupt caches, gzip/payload/HTTP failures, 404 versus 403, transient retry, ambiguous report records, author/model/insertion matching, no metric imputation, full-universe audit preservation and stale-success invalidation.

A real official 11gl XML download/probe succeeds; every one of its eight mapped SG options matches a report record. Its ASU neighbor counts match independent brute-force atom distances. An AlphaFold sample retains 21 context-pending local-pass and 14 local-fail/unassessed options. The end-to-end representative integration uses those two accepted coordinate files and their auxiliary inputs, retains all 760 proteins / 10,799 original site statuses, reuses the real XML cache when available and leaves all review rows pending. Full 1,176-file / 445-report workstation execution is pending. Compact/shared atom indexing stores only context fields; complete ASU atom context is still potentially memory-intensive for very large complexes.

## Primary sources

- Official validation reports: https://www.wwpdb.org/validation/validation-reports
- Official legacy report directory example: https://files.wwpdb.org/pub/pdb/validation_reports/1g/11gl/
- Gemmi biological assemblies/subchain recipes: https://gemmi.readthedocs.io/en/stable/mol.html#assembly
- wwPDB chemical component dictionary (standard component atom names): https://www.wwpdb.org/data/ccd
