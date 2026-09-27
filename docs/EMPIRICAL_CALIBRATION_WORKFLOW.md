# Empirical Calibration Workflow — NEC–FPT

## Step 1 — Inspect all bulk reference datasets

Run:
```r
source("scripts/00_inspect_reference_datasets.R")
```

The inspection must precede DE analysis.

## Step 2 — Start with count-ready FPT calibration

### GSE182638
Two human cell lines:
- MM1R
- MM1S

Per cell line:
- RSL3 vs untreated
- RSL3 vs RSL3+Fer-1

Compute rescue-validated evidence, then aggregate cell lines within study.

### GSE255459
Three human cell lines:
- BT549
- HS578
- SUM159

Per cell line:
- Erastin vs DMSO
- RSL3 vs DMSO

Aggregate inducer evidence within cell line, then cell lines within study.

### GSE319384
Do **not** use for initial weight fitting unless needed.
Prefer as held-out orthogonal FPT validation.

## Step 3 — Reprocess primary NEC data

### GSE108621

Final count-level analysis requires reprocessing raw SRA.

Per gene compute:
1. TSZ vs DMSO
2. TSZ vs TNF
3. TSZ vs TSZ+Nec-1s

Then:

```text
support =
 min(
   max(e_TSZ_DMSO,0),
   max(e_TSZ_TNF,0),
   max(e_Nec1s_reversal,0)
 )

contradiction =
 max(
   max(-e_TSZ_DMSO,0),
   max(-e_TSZ_TNF,0),
   max(-e_Nec1s_reversal,0)
 )

e_GSE108621 = support-contradiction
```

This requires a feature to be:
- induced in necroptosis;
- stronger than TNF inflammation alone;
- reversible by Nec-1s.

### Replication
Next:
- GSE172027
- GSE154230

## Step 4 — Orthogonal NEC validation

Hold:
- GSE134234
- GSE268650

These use direct RIPK3 activation and reduce dependence on TSZ-specific transcription.

## Step 5 — Per-contrast effect engine

For count data:
- DESeq2 MLE effect + SE
- adaptive shrinkage with ashr
- local false sign rate

No hard p-value cutoff defines evidence.

```text
Magnitude = percentile_rank(|posterior beta|)
Confidence = 1-lfsr
e = sign(d*beta)*Magnitude*Confidence
```

## Step 6 — Anti-circularity

If a gene is directly manipulated:
- exclude that experiment from its own E.

Examples:
- RIPK3 in engineered RIPK3 activation studies;
- PSAP in PSAP-KO neurons.

## Step 7 — Cross-study E

Each independent study contributes at most one signed evidence value per gene.

```text
Positive      = mean(max(e_study,0))
Contradiction = mean(max(-e_study,0))
E             = max(0,Positive-Contradiction)
```

## Step 8 — Specificity S

For NEC genes:
- evaluate same-direction behavior in FPT and generic stress.

For FPT genes:
- evaluate same-direction behavior in NEC and generic stress.

```text
C = max(competing_same_direction,stress_same_direction)
S = max(0,(E_target-C)/(E_target+C+epsilon))
```

## Step 9 — Final weights

```text
W_MPS  = M*E
W_NFSI = M*E*S
```

## Step 10 — ESR

Build separately:
- rescue-validated NEC empirical signature
- rescue/cross-inducer FPT empirical signature

Primary scoring should be single-sample and rank-based.

## Step 11 — Required audit outputs

Every study analysis writes:
- metadata.csv
- QC summary
- filtered gene universe
- per-contrast effects
- per-contrast e
- within-study aggregate e
- sessionInfo
- analysis decision notes

Cross-study stage writes:
- E
- same-direction competitor/stress evidence
- S
- final weights
- replication depth.
