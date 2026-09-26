# Protocol Blueprint

## 1. Objective

Develop a reproducible transcriptome-based framework that distinguishes **disulfidptosis (DPT)** from **ferroptosis (FPT)** while separately measuring:

- mechanistic permissiveness,
- mechanistic completeness,
- empirical resemblance to experimentally confirmed states,
- and DPT-vs-FPT specificity.

The framework is deliberately disease-agnostic.

---

## 2. Evidence model

Each candidate gene is assigned:

- program: DPT or FPT
- mechanistic stage
- mechanistic role
- expected RNA direction
- evidence tier
- mechanistic prior `M`
- MPS eligibility
- DSI eligibility
- empirical reproducibility `E`
- specificity `S`
- final MPS weight
- final DSI weight

### Evidence tiers

| Tier | Operational definition | Mechanistic prior M |
|---|---|---:|
| A+ | Direct perturbation plus mechanistic phenotype and/or rescue | 1.00 |
| A | Direct functional perturbation with relevant phenotype | 0.80 |
| B | Strong screen/proteomics/pathway evidence with limited gene-specific functional validation | 0.50 |
| C | Context/pathway evidence without strong direct functional validation | 0.25 |

These are **priors**, not final weights.

---

## 3. DPT architecture

### D1 — Cystine-loading permissiveness

Core:
- SLC7A11
- SLC3A2

Interpretation: capacity for cystine loading. This stage is critical for DPT permissiveness but is not DPT-specific by itself.

### D2 — Reducing-capacity vulnerability

Core:
- G6PD — protective, inverse direction
- TXNRD1 — protective, inverse direction

Extended:
- PGD
- TKT
- TALDO1
- SLC2A1
- SLC2A3

Interpretation: lower glucose/PPP/thioredoxin reducing support increases susceptibility to disulfide stress.

### D3 — WRC/actin execution permissiveness

Core:
- NCKAP1
- WASF2
- CYFIP1
- ABI2
- BRK1

Mechanistic annotation:
- RAC1 activity is mechanistically important, but RAC1 mRNA is not treated as a direct activity readout.

### D4 — Empirical DPT State Resemblance

This is learned from experimentally confirmed DPT perturbation transcriptomes and is not defined by a hand-picked literature list.

### DPT target-availability annotation

Examples:
- ACTB
- MYH9
- MYH10
- FLNA
- FLNB
- TLN1
- ACTN4
- IQGAP1
- MYL6
- DSTN
- CAPZB
- CD2AP
- PDLIM1

These are **not** treated as direct RNA activity markers because the defining DPT evidence is primarily abnormal protein disulfide bonding and actin-cytoskeleton collapse rather than transcriptional induction.

---

## 4. FPT architecture

### F1 — Iron availability

Core:
- NCOA4
- TFRC

Extended protective/context:
- FTH1
- FTL
- SLC40A1

### F2 — PUFA-phospholipid susceptibility

Core:
- ACSL4
- LPCAT3

### F3 — Lipid-peroxidation machinery

Core:
- POR

Extended:
- ALOX15

### F4 — Anti-ferroptotic defense failure

Core protective genes scored inversely:
- GPX4
- AIFM2 / FSP1
- DHODH
- GCH1

Extended:
- DHFR

---

## 5. Shared biology

The following can contribute strongly to biological context but should not directly dominate DPT-vs-FPT specificity:

- SLC7A11 / SLC3A2
- GSH-related machinery
- NRF2 response
- generic ROS response
- broad NADPH/redox programs

Shared components remain usable for permissiveness and interaction terms.

---

## 6. Reference perturbation library

Reference data are divided into:

1. experimentally supported DPT perturbations
2. experimentally supported FPT perturbations
3. generic/competing stress controls

The seed library is stored in `config/reference_datasets.csv`.

### DPT factorial design

For TXNRD1 loss × glucose deprivation, avoid a naïve final-state contrast because it mixes genotype effects, glucose-starvation effects, and DPT.

Use the interaction:

```text
(KO_GluMinus - KO_Regular) - (EV_GluMinus - EV_Regular)
```

This aims to isolate the response specific to reducing-defense failure under glucose deprivation.

### FPT rescue design

For an inducer + Ferrostatin-1 dataset:

- induction contrast: inducer vs control
- rescue contrast: inducer vs inducer + Ferrostatin-1

High-confidence empirical FPT behavior should be induced by the FPT trigger and reversed by rescue.

---

## 7. Anti-circularity rule

A directly manipulated gene cannot use the same experiment to validate its own RNA behavior.

Examples:
- TXNRD1 KO does not validate TXNRD1 as a transcriptomic DPT marker.
- GPX4 KO does not validate GPX4 as a transcriptomic FPT marker.
- SLC7A11 KO does not validate SLC7A11 as a transcriptomic marker.

This rule is mandatory.

---

## 8. Empirical per-study evidence

For gene `g` in study `j`:

1. estimate a shrunken effect size `beta(g,j)`
2. estimate sign confidence, preferably local false-sign rate (`lfsr`)
3. rank `|beta|` within the same study

```text
Magnitude(g,j)  = percentile_rank(|beta(g,j)|)
Confidence(g,j) = 1 - lfsr(g,j)
```

Let `d(g)` be the expected direction:
- +1 = higher RNA supports the state
- -1 = protective gene; lower RNA supports the state
- 0 = do not score RNA direction directly

```text
e(g,j) =
  sign[d(g) × beta(g,j)]
  × Magnitude(g,j)
  × Confidence(g,j)
```

Range: -1 to +1.

Positive values support the expected biology.
Negative values are contradictions and must not be silently discarded.

If several contrasts are from one study, aggregate those contrasts first so one study does not dominate replication.

---

## 9. Empirical reproducibility score E

Across independent studies:

```text
Positive(g)      = mean(max(e(g,j), 0))
Contradiction(g) = mean(max(-e(g,j), 0))

E(g) = max(0, Positive(g) - Contradiction(g))
```

Range: 0 to 1.

Independent-study count is reported separately as replication depth.

---

## 10. Specificity score S

For a DPT gene, compare target-state evidence against **same-direction** behavior in FPT and generic stress datasets.

```text
C(g) = max(
  FPT_same_direction(g),
  Stress_same_direction(g)
)
```

```text
S(g) = max(
  0,
  [E_target(g) - C(g)] /
  [E_target(g) + C(g) + epsilon]
)
```

Range: 0 to 1.

The same logic is mirrored for FPT.

Opposite-direction behavior in a competing program is not treated as mimicry.

---

## 11. Final gene weights

Permissiveness:

```text
W_MPS(g) = M(g) × E(g)
```

Specificity:

```text
W_DSI(g) = M(g) × E(g) × S(g)
```

A gene can therefore remain mechanistically important while receiving little RNA diagnostic weight.

---

## 12. Single-sample scoring

Within each sample, rank genes independently:

```text
r(g,s) in [0,1]
```

For promoter genes:

```text
q(g,s) = r(g,s)
```

For protective genes:

```text
q(g,s) = 1 - r(g,s)
```

Thus `q` near 1 always means stronger support for the corresponding death-state mechanism.

---

## 13. Stage scores

```text
Stage_k(s) =
  sum[W(g) × q(g,s)] /
  sum[W(g)]
```

Range: 0 to 1.

---

## 14. DPT aggregation

### Initiation gate

```text
DPT_init = sqrt(D1 × D2)
```

This prevents high cystine-loading evidence from being interpreted as strong DPT permissiveness when reducing-capacity vulnerability is absent.

### Molecular permissiveness

```text
DPT-MPS = (DPT_init + D3) / 2
```

### Mechanistic completeness

```text
DPT-MCI = (D1 × D2 × D3)^(1/3)
```

MCI is deliberately weakest-link sensitive.

---

## 15. FPT aggregation

```text
FPT_init = (F1 × F2 × F3)^(1/3)

FPT-MPS =
  (FPT_init + F4) / 2

FPT-MCI =
  (F1 × F2 × F3 × F4)^(1/4)
```

---

## 16. Empirical State Resemblance

Keep empirical transcriptomic resemblance separate from mechanistic permissiveness:

- DPT-ESR
- FPT-ESR

Primary implementation should be rank-based and single-sample. Alternative methods such as ssGSEA/GSVA may be used for sensitivity analysis.

---

## 17. Disulfidptosis Specificity Index

Build `D_specific` and `F_specific` from high-specificity evidence only.

Shared redox components are excluded or strongly down-weighted.

```text
DSI =
  (D_specific - F_specific) /
  (D_specific + F_specific + epsilon)
```

Range: -1 to +1.

- positive = DPT-skewed
- negative = FPT-skewed
- near zero = no mechanism-specific dominance

DSI is **not a probability**.

---

## 18. Classification rules

Each program is considered supported only when its calibrated criteria for:
- MPS
- MCI
- ESR

are satisfied.

| DPT support | FPT support | Final class |
|---|---|---|
| yes | no | DPT-dominant |
| no | yes | FPT-dominant |
| yes | yes | Mixed |
| no | no | Indeterminate / Neither |

For Mixed samples, DSI may annotate:
- Mixed — DPT-skewed
- Mixed — balanced
- Mixed — FPT-skewed

Partial states are explicitly reportable:
- partial DPT permissiveness
- incomplete FPT mechanism

---

## 19. Threshold calibration

Thresholds are not hard-coded.

Recommended:
- ROC-based candidate thresholds
- bootstrap threshold stability
- leave-one-study-out validation

The held-out unit must be an **entire study**, not random samples from the same experiment.

If a threshold is unstable, retain the score as continuous evidence and mark classification as insufficiently calibrated.

---

## 20. Required output

For every scored sample:

- D1, D2, D3
- DPT-MPS
- DPT-MCI
- DPT-ESR
- F1, F2, F3, F4
- FPT-MPS
- FPT-MCI
- FPT-ESR
- D_specific
- F_specific
- DSI
- final class
- reference-depth/confidence annotation

---

## 21. Freeze policy

Before any future application:

- freeze gene dictionaries
- freeze stage definitions
- freeze equations
- freeze calibration procedure
- version every later change

No future application dataset may be used to redesign the classifier after its results are seen.
