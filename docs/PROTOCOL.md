# Protocol Blueprint — NEC–FPT Framework

## 1. Objective

Develop a reproducible transcriptome-based framework that distinguishes **necroptosis (NEC)** from **ferroptosis (FPT)** while separately measuring:

- mechanistic permissiveness,
- mechanistic completeness,
- empirical resemblance to confirmed perturbation states,
- and mechanism-specific dominance.

The framework is disease-agnostic.

---

## 2. Evidence model

Each candidate gene is assigned:

- program
- mechanistic stage
- mechanistic role
- expected RNA direction
- evidence tier
- mechanistic prior `M`
- MPS eligibility
- NFSI eligibility
- empirical reproducibility `E`
- specificity `S`
- final MPS weight
- final NFSI weight

### Evidence tiers

| Tier | Operational definition | M |
|---|---|---:|
| A+ | Direct perturbation plus mechanistic phenotype and/or rescue | 1.00 |
| A | Direct functional perturbation with relevant phenotype | 0.80 |
| B | Strong pathway/screen evidence with limited gene-specific validation | 0.50 |
| C | Context association without strong direct validation | 0.25 |

These are priors, not final weights.

---

## 3. Necroptosis architecture

### N1 — Necrosome competence

Core:
- **RIPK3**
- **RIPK1**

Interpretation:
- RIPK3 is the central necroptotic kinase.
- RIPK1 is important in canonical TNFR1-driven necroptosis but is not universally required because RIPK3 can also be activated through other RHIM-containing routes.

RNA interpretation:
higher expression may increase molecular competence, but phosphorylation/complex formation is required for execution.

### N2 — MLKL execution competence

Core:
- **MLKL**

Interpretation:
RIPK3-mediated MLKL activation, oligomerization and membrane localization form the terminal execution module.

RNA interpretation:
MLKL expression indicates execution capacity; it does **not** demonstrate pMLKL or MLKL oligomerization.

### Trigger-context / route annotations

Extended:
- ZBP1
- TICAM1
- TNFRSF1A
- CYLD

These provide upstream route context but do not define necroptosis by themselves.

### Checkpoint annotations not directly scored from bulk RNA

- CASP8
- FADD
- CFLAR

Reason:
their necroptosis-related effects depend heavily on catalytic state, complex composition and/or isoform balance. Bulk transcript direction is not safely interpretable as execution activity.

### NEC empirical state

`NEC-ESR` is learned from experimentally confirmed necroptosis perturbation transcriptomes, especially rescue-validated and orthogonal induction designs.

---

## 4. Ferroptosis architecture

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

### F4 — Anti-ferroptotic defense
Core protective genes scored inversely:
- SLC7A11
- GPX4
- AIFM2 / FSP1
- DHODH
- GCH1

Extended:
- SLC3A2
- DHFR

---

## 5. Mechanistic boundary

The primary distinction is:

```text
NEC:
RIPK3 -> MLKL -> membrane permeabilization

FPT:
iron-dependent PUFA-phospholipid peroxidation
+ failure of lipid-antioxidant defense
```

Generic inflammatory or oxidative stress is not sufficient to classify either state.

---

## 6. Reference perturbation library

Reference data are divided into:

1. rescue-validated NEC calibration studies
2. orthogonal NEC validation studies
3. rescue/cross-inducer FPT calibration studies
4. orthogonal/neural FPT validation studies
5. generic stress comparators

See `config/reference_datasets.csv`.

### Primary NEC calibration: GSE108621

Use three aligned comparisons:
- TSZ vs DMSO
- TSZ vs TNF
- TSZ vs TSZ+Nec-1s

A strong NEC empirical feature should:
1. appear under TSZ;
2. exceed TNF-only inflammatory behavior;
3. be reversed by Nec-1s.

### NEC replication

GSE172027:
- TSZ vs DMSO
- TSZ vs TSZ+Nec-1s
in human astrocytes.

GSE154230:
- microglia and astrocytes;
- TSZ and an additional RIPK1-activating route;
- Nec-1s rescue.

### Orthogonal NEC validation

GSE134234:
direct chemical dimerization of engineered RIPK3 with TNF-only inflammatory control.

GSE268650:
optogenetic RIPK3 oligomerization with dark/light and RIPK3-inhibitor controls.

GSE287439:
single-cell human iPSC neuron/astrocyte/microglia tri-culture for cell-type-specific validation.

### FPT calibration

Primary rescue:
- GSE182638: two human lines, RSL3 ± Ferrostatin-1.

Independent rescue:
- GSE247883: A549, RSL3 ± Ferrostatin-1.

Cross-inducer:
- GSE255459: Erastin and RSL3 across three human lines.

Orthogonal held-out validation:
- GSE319384: Erastin, RSL3 and Ferroptocide in isogenic PANC-1 backgrounds.

Neural validation:
- GSE287284
- GSE152988

---

## 7. Anti-circularity

A directly manipulated gene cannot use the same experiment to validate its own transcriptomic behavior.

Examples:
- engineered RIPK3 activation datasets do not validate RIPK3 RNA as an empirical marker of itself;
- a GPX4 knockout experiment cannot validate GPX4 transcript behavior;
- PSAP-KO cannot validate PSAP itself.

---

## 8. Per-contrast empirical evidence

For gene `g` in contrast `j`:

```text
Magnitude(g,j)  = percentile_rank(|beta(g,j)|)
Confidence(g,j) = 1 - lfsr(g,j)

e(g,j) =
  sign[d(g) * beta(g,j)]
  * Magnitude(g,j)
  * Confidence(g,j)
```

Range: [-1,+1].

Negative evidence is retained as contradiction.

---

## 9. Rescue-validated evidence

For an induction/rescue pair:

```text
support =
  min(max(e_induction,0),
      max(e_reversal,0))

contradiction =
  max(max(-e_induction,0),
      max(-e_reversal,0))

e_rescue = support - contradiction
```

### Three-way NEC control in GSE108621

For:
- TSZ vs DMSO
- TSZ vs TNF
- TSZ vs TSZ+Nec-1s

use:

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

e_rvNEC = support - contradiction
```

This strongly penalizes generic TNF inflammatory genes.

---

## 10. Across-study empirical reproducibility E

Each independent study contributes at most one signed value per gene.

```text
Positive(g)      = mean(max(e_study,0))
Contradiction(g) = mean(max(-e_study,0))

E(g) = max(0, Positive(g) - Contradiction(g))
```

Study count is reported separately as replication depth.

---

## 11. Specificity S

For a NEC gene, evaluate the **same NEC direction** in FPT and generic-stress datasets.

For an FPT gene, evaluate the **same FPT direction** in NEC and generic-stress datasets.

```text
C(g) = max(
  competing_program_same_direction,
  generic_stress_same_direction
)

S(g) =
  max(
    0,
    (E_target(g)-C(g)) /
    (E_target(g)+C(g)+epsilon)
  )
```

Opposite-direction behavior in the competing program is not mimicry.

---

## 12. Final weights

```text
W_MPS(g)  = M(g) * E(g)
W_NFSI(g) = M(g) * E(g) * S(g)
```

---

## 13. Single-sample scoring

Within each sample:

```text
r(g,s) in [0,1]
```

Promoter:
```text
q(g,s)=r(g,s)
```

Protective:
```text
q(g,s)=1-r(g,s)
```

Stage:
```text
Stage_k(s)=sum(W_g*q_gs)/sum(W_g)
```

---

## 14. NEC aggregation

```text
N1 = Necrosome competence
N2 = MLKL execution competence

NEC-MPS = (N1 + N2) / 2
NEC-MCI = sqrt(N1 * N2)
```

MPS captures overall permissiveness.

MCI is weakest-link sensitive: high RIPK1/RIPK3 competence without MLKL execution competence, or vice versa, cannot produce a high complete-mechanism score.

`NEC-ESR` remains independent.

---

## 15. FPT aggregation

```text
FPT_init = (F1 * F2 * F3)^(1/3)
FPT-MPS  = (FPT_init + F4) / 2
FPT-MCI  = (F1 * F2 * F3 * F4)^(1/4)
```

`FPT-ESR` remains independent.

---

## 16. Necroptosis–Ferroptosis Specificity Index

```text
NFSI =
  (NEC_specific - FPT_specific) /
  (NEC_specific + FPT_specific + epsilon)
```

Range [-1,+1]:
- positive = NEC-skewed
- negative = FPT-skewed
- near zero = no clear mechanism-specific dominance

NFSI is not a probability.

---

## 17. Final classes

| NEC support | FPT support | Class |
|---|---|---|
| yes | no | NEC-dominant |
| no | yes | FPT-dominant |
| yes | yes | Mixed |
| no | no | Indeterminate / Neither |

Support requires calibrated MPS, MCI and ESR criteria.

---

## 18. Calibration

No manual final thresholds.

Use:
- reference-positive/negative distributions;
- bootstrap threshold stability;
- leave-one-study-out validation;
- whole-study holdout rather than random sample split.

Unstable thresholds remain continuous evidence rather than forced classes.

---

## 19. Freeze policy

Before any external application:
- freeze dictionaries;
- freeze equations;
- freeze empirical signatures;
- freeze thresholds;
- freeze software version.

No external application dataset may be used to redesign the classifier after results are seen.
