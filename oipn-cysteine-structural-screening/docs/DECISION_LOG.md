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
