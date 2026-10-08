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
