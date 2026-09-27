# Ferroptosis-Only Protocol — Draft v0.1

## 1. Goal

Develop a reproducible transcriptomic model of ferroptosis that is mechanistically anchored and empirically calibrated.

The model is **not** intended to prove lipid peroxidation from RNA alone. Transcriptomic scores represent molecular permissiveness, completeness and resemblance to experimentally confirmed ferroptotic states.

## 2. Initial mechanistic modules

### F1 — Iron availability
Candidate biology:
- ferritinophagy
- iron uptake
- iron storage/export

### F2 — PUFA-phospholipid susceptibility
Candidate biology:
- PUFA activation
- PUFA incorporation into membrane phospholipids

### F3 — Lipid-peroxidation machinery
Candidate biology:
- enzymatic and non-enzymatic drivers of phospholipid oxidation

### F4 — Anti-ferroptotic defense
Candidate biology:
- System Xc− / GSH / GPX4
- FSP1 / CoQ
- DHODH / mitochondrial CoQ
- GCH1 / BH4

### F5 — Regulatory/context layer
Genes with clear ferroptosis biology but insufficient specificity or unsafe direct RNA interpretation are retained as context rather than forced into the core score.

## 3. Evidence classes

- **CORE** — direct mechanistic importance and suitable for RNA scoring
- **EXTENDED** — validated regulator/context, but not universally required
- **MECHANISTIC** — biologically important but unsafe for direct transcript-direction scoring
- **EMPIRICAL_ONLY** — reproducible ferroptotic transcriptomic marker without sufficient causal status for the mechanistic core
- **EXCLUDED** — insufficient or contradictory evidence

## 4. Empirical calibration

Reference perturbation studies should prioritize:
- canonical ferroptosis induction
- rescue by Ferrostatin-1/Liproxstatin-1 where available
- multiple inducers
- multiple cell lines
- neuronal models

A directly manipulated gene cannot validate itself in the same experiment.

## 5. Primary outputs

### FPT-MPS
Molecular permissiveness score.

### FPT-MCI
Mechanistic completeness score.

### FPT-ESR
Empirical state resemblance score.

These remain separate rather than being collapsed into one opaque score.

## 6. Negative controls

To avoid building an oxidative-stress detector, validation should include:
- H2O2 oxidative stress
- non-ferroptotic cytotoxic stress
- ideally apoptosis/necroptosis perturbation references in later phases

## 7. Freeze policy

The gene dictionary and equations must be frozen before external biological application datasets are analyzed.
