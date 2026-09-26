# Empirical Calibration Workflow

This document specifies the executable order for Milestones 1–4.

## Step 1 — Ingest count-ready studies first

Order:

1. GSE255459 — FPT, three human cell lines, Erastin/RSL3
2. GSE131444 — FPT, mouse MEF, Erastin
3. GSE317656 — FPT, two human neuroblastoma lines, Erastin
4. GSE104664 — generic oxidative stress, H2O2 time course

These four datasets can establish and test the count-based analysis engine before raw-SRA reprocessing is introduced.

## Step 2 — Per-study QC

For every study:

- reconstruct sample metadata from GEO;
- verify integer count matrix;
- confirm column/sample identity;
- report library sizes;
- report number of detected genes;
- filter low-count genes using the locked rule in `config/analysis_parameters.csv`;
- perform VST PCA and sample correlation;
- record, but do not silently remove, any potential outlier.

No raw expression matrices are merged across studies.

## Step 3 — Differential effects

Use DESeq2.

For each pre-specified contrast, obtain:
- MLE log2 fold change;
- standard error.

Then apply adaptive shrinkage using `ashr`:
- posterior mean effect (`beta`);
- local false sign rate (`lfsr`).

No hard p-value/FDR cutoff is used to define empirical evidence.

## Step 4 — Convert each contrast into evidence e

For target direction d:

```text
Magnitude = percentile_rank(|posterior beta|)
Confidence = 1 - lfsr

e = sign(d × beta) × Magnitude × Confidence
```

Range: [-1, +1].

## Step 5 — Within-study aggregation

### GSE255459

For each cell line:

```text
e_cell = median(e_Erastin, e_RSL3)
```

Then:

```text
e_study = median(e_BT549, e_HS578, e_SUM159)
```

This ensures the study contributes one evidence unit despite six contrasts.

### GSE317656

```text
e_study = median(e_SK-N-AS, e_KELLY)
```

### GSE104664

```text
e_study = median(e_H2O2_16h, e_H2O2_36h)
```

### GSE247883 (after raw reprocessing)

Define:
- induction = RSL3 - DMSO
- reversal = RSL3 - (RSL3 + Fer-1)

Both are evaluated using the same FPT target direction.

```text
support = min(max(e_induction,0), max(e_reversal,0))
contradiction = max(max(-e_induction,0), max(-e_reversal,0))
e_rescue = support - contradiction
```

Thus a transcript only receives strong rescue-validated evidence when the RSL3 change is also reversed by Fer-1.

### GSE282334 (after raw reprocessing)

Use the interaction coefficient:

```text
(KO_GluMinus - KO_Regular) - (EV_GluMinus - EV_Regular)
```

TXNRD1 is excluded from self-validation in this study.

## Step 6 — Cross-study E

After each study has exactly one signed evidence value per gene:

```text
Positive       = mean(max(e_study, 0))
Contradiction  = mean(max(-e_study, 0))
E              = max(0, Positive - Contradiction)
```

Independent-study count is stored separately.

## Step 7 — Specificity S

For a DPT target gene, competitor FPT and generic-stress transcriptomes are re-evaluated using the **DPT expected direction**.

For an FPT target gene, DPT and stress datasets are re-evaluated using the **FPT expected direction**.

This is essential: opposite-direction behavior in a competing death program is discriminatory evidence, not mimicry.

```text
C = max(competitor_same_direction, generic_stress_same_direction)

S = max(0, (E_target - C) / (E_target + C + epsilon))
```

## Step 8 — Final gene weights

```text
W_MPS = M × E
W_DSI = M × E × S
```

## Step 9 — Required audit files

Each run must write:
- sample metadata;
- QC summary;
- filtered gene universe;
- per-contrast effect table;
- per-study signed evidence;
- cross-study E;
- same-direction competitor/stress evidence;
- S;
- final gene weights;
- package/session information.
