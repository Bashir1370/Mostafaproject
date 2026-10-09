# Decision log — 2026-10-08

| ID | Decision | Status |
|---|---|---|
| D01 | Mouse replaces rat; mouse sequence/structure is primary | Approved |
| D02 | Oxaliplatin vs matched vehicle in DRG; OIPN-specific initial scope | Approved |
| D03 | Discovery uses significant protein-coding DEGs, both directions, BH FDR <0.05, no FC cutoff | Approved; supersedes the earlier all-expressed-gene proposal |
| D04 | Ten scientific steps: 01–02 input, 03–09 structural/evaluation, 10 interpretation | Approved |
| D05 | Unit of scoring is Cys; protein summary is maximum site score | Approved |
| D06 | Two equal-weight percentile features: SG accessibility and -pKa | Frozen initial heuristic; biological performance untested |
| D07 | Cys–Cys distance supplemental; Cys count, helices, FDR, FC and function not score terms | Approved |
| D08 | Target predicted Cys pLDDT >=90; neighbors >=70 within 6 Å; stricter sensitivity analysis | Approved pragmatic QC |
| D09 | Incompatible chemical states and missing data separated; no zero susceptibility imputation | Approved |
| D10 | Empirical top-quartile shortlist with boundary ties retained; all rankings exported | Approved prioritization rule |
| D11 | Independent evaluation cannot be replaced by enrichment or rank changes alone | Approved |
| D12 | Ubuntu/Bash workflow; independent project configuration/results | Approved |
| D13 | New folder oipn-cysteine-structural-screening; no modification of existing scientific workflows | Approved |
| D14 | GSE286387 primary candidate; GSE125002 complementary candidate | Pending current Step 01 audit |

## Pending implementation choices, not biological result claims

Exact Ubuntu/R/Bioconductor/Python versions; dataset/sample accession lock; model terms dictated by audited design; reference releases; experimental local-validation implementation and deterministic structure tie selection; SASA implementation/radii and PROPKA version; benchmark supplement availability and experiment-derived labels; enrichment term collection. These are resolved and recorded before the corresponding analysis. No sample GSM IDs are invented in advance.

## Relation to the motivating article

Reuse its sequence–structure integration, systematic filtering, transparent ranking and downstream interpretation architecture. Do not copy histidine helix weights or the tumor pH window. Our working feature directions are hypotheses for prioritization, with a separate evaluation gate. We are not claiming a specific S-glutathionylation predictor or measured oxidation in OIPN.

## Amendment policy

Version every parameter/method change and state why. Preserve v0.1.0 results if a later model is trained. Do not change the initial rule retrospectively to obtain a desired named candidate.

## 2026-10-08 — Step 01 input implementation

Official SOFT and submitter-designated raw counts retrieved and audited. Candidate GSE286387 remains provisional until R sample QC and animal/pooling review. Retain one vector for duplicate stable IDs only if every count and biotype matches; preserve annotation alternatives, never sum copies. Actual 167 duplicate IDs / 276 redundant rows satisfy this rule. No scientific threshold changed. A protocol-independent input integrity script uses Python standard library; QC and DE use R next. Sample identities are explicit; batch/age/animal IDs remain NA rather than inferred.

## 2026-10-08 — Step 01b sample QC implementation

User workstation reports R 4.3.3 / Bioconductor 3.18 / DESeq2 1.42.1. QC uses an intercept-only, blind DESeq2 VST before any Wald/LRT testing. PCA uses top 500 variable retained genes without feature scaling; Euclidean distances and Pearson correlations use all retained VST genes. No sample is automatically excluded. Count filter and discovery thresholds are unchanged. Exact source checksums/identity mapping and every processed count are checked against original deposited counts. Full seven-test suite and actual-data run pass in R 4.3.3 / DESeq2 1.42.0; user reproduction with 1.42.1 pending. Metadata remains held for animal/pooling/design review.

## Step 01 acceptance and Step 02 implementation — 2026-10-08

The user's actual Ubuntu QC snapshot is committed and reviewed. Study methods and the reported ten-mouse RNA-seq design support five independent units per group; library pooling does not establish cross-animal tissue pooling. Original animal IDs, batch, collection age and exact final-dose interval remain documented unknowns. Retain all ten samples, including DRG_10mg_rep1; no independent technical evidence supports its exclusion. Freeze `~ condition`, control reference, oxaliplatin numerator, with accepted-manifest/count checksums. See DATASET_DESIGN_AUDIT.md for evidence and the distinction between reported facts and our design interpretation.

Step 02 uses the existing registered count filter, BH across valid tested genes of all biotypes, independentFiltering=FALSE, default Cook's handling, and coding padj<0.05 selection without FC/baseMean cutoffs. Real-input implementation validation (R 4.3.3 / DESeq2 1.42.0) and all eleven tests passed. The user's 1.42.1 DE reproduction remains pending. No scoring thresholds changed; protocol remains 0.1.0. The earlier provisional/held entries are historical stages superseded by this acceptance.

## Step 02 workstation report — 2026-10-08

The user supplied de_report.md from Ubuntu with R 4.3.3 / DESeq2 1.42.1 and STEP02_COMPLETE. All reported totals match the separate 1.42.0 implementation run: 18,538 tested, 18,503 valid p-values, 893 significant across biotypes, 854 coding (477 up / 377 down). Preserve the supplied report verbatim under docs/audits/ubuntu_DE_2026-10-08. Actual gene identities, checksum records and environment files must be uploaded from the workstation and inspected before accepting the Step 03 input. No local 1.42.0 result is relabeled as a 1.42.1 file. Parameters and design are unchanged.

## Step 02 snapshot acceptance — 2026-10-08

Reviewed actual workstation files from d80a26859ad54577ebf4021e87ef4f4a0cb8cac0. All eight SHA256 checks, frozen input/configuration records, ten-sample identities and design matrix pass. Discovery CSV contains 854 unique significant coding IDs (477 up / 377 down); all IDs match the separate validated 1.42.0 discovery run. Record review limits and checksum-bound Step 03 input. Full workstation BH was not independently recomputed because all-gene p-values are outside this selected snapshot. Protein mapping has not started; scientific thresholds and protocol version are unchanged.

## Step 03 mapping implementation — 2026-10-08

Operationalize the frozen ambiguous-mapping hold policy with exact current Ensembl GeneId cross-references in UniProtKB, taxonomy 10090, all reviewed/unreviewed entries retained as candidates, and acceptance only of a unique accession with compatible displayed-sequence association and valid nonfragment sequence. Numeric ID version suffixes are recorded then removed for stable-ID comparison; no claim of unchanged release-76 sequence follows. Preserve both explicit displayed-isoform and entry-level unspecified associations. Do not choose among multiple entries by review status, length, function or eventual structural score. Hold missing current links for separately logged historical-ID rescue rather than using symbols. Record current reference release/access dates, canonical/precursor limitations and many-gene/one-accession cases.

Real accepted-input validation: UniProt 2026_03, 285 passing genes / 283 unique proteins; 562 multi-entry holds, six no-current-link holds, one fragment unassessable. All nineteen tests pass. This first pass leaves a substantial ambiguous cohort and is not the complete candidate proteome. Workstation reproduction and identity-resolution review precede Step 04. These implementation choices are fixed before structural candidate scoring; protocol 0.1.0 and existing numeric thresholds are unchanged. Existing scientific workflows and frozen DE inputs/configuration are preserved.

## Step 03b representative resolution and coverage objective — 2026-10-08

User requires mapping most genes before cysteine inventory. Do not turn this into a forced percentage target or accept unsupported identities. Add evidence-based resolution version 1.0.0 within scientific Step 03: explicit Ensembl gene links in the EMBL-EBI gene-centric service, one representative, complete original candidate accounting, representative among original exact-gene mouse UniProt candidates, and original canonical/sequence checks. Review status and sequence length are not selection rules. Preserve first-pass evidence and statuses; output a combined mapping in a separate directory. Six no-current-link genes and fragment/remaining representative conflicts remain documented for separate review.

Real-input validation resolves 504 additional genes: 789 / 854 (92.4%), 787 unique proteins, 64 held plus one unassessable fragment. All twenty-six tests pass. Gene-centric service release headers were absent; document raw relationship response hashes and retrieval times, not a fictitious 2026_03 relationship release. Sequences remain on the original UniProt release. Workstation first-pass matching totals are reported, but full workstation mapping files and Step 03b reproduction await review. Existing thresholds, frozen DE inputs and original Step 03 script are unchanged; protocol 0.1.0 remains the study protocol.

## 2026-10-09 — Mapping snapshot review and Step 04

- Reviewed both actual Ubuntu mapping snapshots at source commit 3bdf58061337550aeda2f180e44fc710f810170b. All 20 snapshot hashes, upstream bindings, exact mapping/FASTA/candidate/evidence tables and response body hashes pass. Accept only the 789 mapped genes / 787 accessions; the remaining 65 cases stay open with explicit reasons.
- Freeze Step 04 inputs in a new binding without changing historical analysis_parameters.json, Step 03 input or previous scripts/results. Use the committed reviewed snapshot offline.
- Count only canonical C using 1-based full-sequence numbering, preserve precursor regions and U/O, deduplicate by accession/site, audit zero-C exclusions and all upstream unresolved genes. Counts/percentages remain descriptive.
- Validation execution yields 760 C-positive proteins, 27 zero-C proteins and 10,799 sites. User Ubuntu reproduction is pending; no structural interpretation is inferred.

## 2026-10-09 — Track results directly in Git

User explicitly requested that results/ be committed. Remove only the results/ ignore rule; raw data, caches, environments and logs retain their ignore rules. Update execution documentation to require explicit add/commit/push after runs. Preserve existing reviewed snapshots and checksum-bound scientific inputs. This is a version-control change, with no scientific parameter or analysis-script change. The user uploads their actual workstation results; validation files are not relabeled or uploaded as workstation output.

## 2026-10-09 — Step 04 workstation acceptance and Step 05a

- Accept actual Ubuntu Step 04 snapshot at e419a4e98e8b51368a0d5123b08fc128415eee1c after all nine snapshot / 27 upstream hashes and complete sequence/site/audit reproduction checks. Bind 760 C-positive accessions and 10,799 unique sites without modifying older scientific inputs.
- Split implementation of scientific Step 05 into metadata catalogue (05a) and subsequent coordinates/residue/assembly/local-quality work (05b); the study still has exactly ten scientific steps.
- Query both PDBe/SIFTS and AlphaFold for every accession. Retain all offered records; no first/best-resolution/global-confidence selection. Exact mouse/accession/canonical-sequence metadata gates are predeclared before examining candidate identities. Partial and mismatched/alternate-isoform models retain explicit audit reasons.
- Use API-provided model URLs/versions and current field names with conflict-checked legacy fallback. Cache SHA256-bound response bytes including 404s; network errors fail rather than become missing structures. Metadata does not certify Cys/SG coverage or quality.
- Historical structural thresholds/score weights remain unchanged; no coordinate extraction, chemical-state analysis or susceptibility score is run in 05a. Actual Ubuntu reproduction follows implementation validation.

Live implementation validation and final cache replay complete all 1,520 queries: 94 proteins with mouse PDBe candidates, 731 with exact-full-sequence AlphaFold candidates, 94 with both, 29 without metadata-eligible candidates (27 unassessable, two with unresolved offered metadata held). All 42 regression tests pass. No coordinates or local-quality approval are implied; Ubuntu reproduction remains pending.

## 2026-10-09 — Workstation Step 05a HTTP 403 diagnostic revision

User reports HTTP 403 before a final catalogue. Original exception handling lost endpoint context. Add diagnostic implementation revision 1.0.1: fatal HTTP errors retain requested URL/status and a bounded, allowlisted failure_context.json. Preserve HTTPError compatibility, 404 handling, existing caches, all successful catalogue tables and scientific inputs. Do not classify the workstation failure as absent structures or claim that access is repaired. Same initial endpoints return expected responses in validation; actual workstation endpoint/response evidence is pending.

Diagnostic revision validation: all 44 regression tests pass; final cache replay reproduces all four complete catalogue/query/protein-audit tables byte-for-byte. Workstation 403 cause and resolution remain pending.

## 2026-10-09 — Step 05a frozen reference reproduction

After workstation HTTP 403 and SSL EOF failures, add explicit --reference-cache mode with all 1,520 public API responses from the previously completed validation run. Preserve original acquisition times/hashes, verify the full sequence/query universe and archive, and import into a separate namespace. Force offline execution and label source reuse in summary/report/input provenance. Do not disable TLS or infer network errors as absent structures. Historical validation records and all frozen scientific inputs remain unchanged.

## 2026-10-09 — Actual Ubuntu Step 05a catalogue accepted

Review source commit f8a1a4e8c60593101425ebc09ba55c2081393bba: eight complete outputs, all 42 input checksums and four byte-identical reference-reproduced tables pass. Freeze actual workstation snapshot with a separate Step 05b input binding. Acceptance concerns catalogue provenance and coverage only, not structure or site eligibility. Retain 29 no-eligible-metadata proteins (27 unassessable, two held with unresolved records) in downstream audits. Coordinate acquisition and local-quality implementation are still pending.

## 2026-10-09 — Step 05b1 raw acquisition and archival

Implement coordinate acquisition as a separate substep before residue mapping/local quality. Download every metadata-eligible option, deduplicating full deposited PDB entries and exact AlphaFold file URLs. Preserve ambiguous/excluded metadata and auxiliary gaps. Add four-kind preflight, verified resumable caches, transport diagnostics and a tar.gz with full member read-back checks; archive existence does not imply complete acquisition or structural acceptance. Raw coordinates remain local under ignored data/; compact manifests are tracked. No frozen scientific scoring/input files or unrelated subprojects are changed.

## 2026-10-09 — Diagnose incomplete acquisition and add evidence-preserving retry

Review actual 2,638-row manifest and 517 failures at ed5ae356e169d2ae5be5b657bbce5130de582bc6. Add missing IncompleteRead/remote-disconnect retries, explicit request timeout/attempt controls and bounded HTTP diagnostics. Do not infer HTTP 400 source without its response body; two exact failed PAE URLs work elsewhere. Probe a recorded failed URL per kind before bulk retry; use one worker for workstation continuation. Preserve successful caches and all frozen identity/scientific gates. All 60 regression tests pass; the partial archive receipt is not remotely validated archive content.
