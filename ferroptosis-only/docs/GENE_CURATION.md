# Systematic Ferroptosis Gene Curation

## Purpose

This document defines how genes are admitted into the ferroptosis-only model.

The goal is **not** to create the largest possible ferroptosis gene list. The goal is to construct a mechanistically interpretable universe that can later be empirically calibrated against perturbation transcriptomes.

## Evidence hierarchy

### A+
Direct genetic/pharmacologic perturbation plus mechanistic readout, rescue, lipidomics or orthogonal validation.

### A
Direct functional perturbation demonstrates a ferroptosis-relevant phenotype, but evidence is less universal or more context dependent.

### B
Strong pathway, CRISPR-screen or biochemical evidence with limited generalization.

### C
Association/context only. C-tier genes are not eligible for the initial mechanistic core.

## Status definitions

### CORE
Required or major ferroptosis determinant with strong mechanistic evidence and a plausible direct RNA interpretation.

### EXTENDED
Validated ferroptosis regulator, but:
- not universally required;
- strongly context dependent;
- or less specific at the RNA level.

### MECHANISTIC
Important to ferroptosis biology, but transcript abundance is not a safe proxy for functional state.

### EMPIRICAL_ONLY
Useful transcriptomic response marker but not a causal core regulator.

### PENDING
Biologically plausible or supported in restricted contexts; must pass empirical calibration before scoring.

## Mechanistic modules

### F1 — Iron availability
Core question:
> Is sufficient redox-active iron available to sustain ferroptotic phospholipid oxidation?

Priority:
- NCOA4
- TFRC

Extended/context:
- IREB2
- FTH1/FTL
- SLC40A1
- SLC11A2

Autophagy genes are retained as mechanistic context, not direct ferroptosis RNA markers.

### F2 — PUFA-phospholipid susceptibility
Core question:
> Does the membrane lipid composition provide oxidizable ferroptotic substrates?

Core:
- ACSL4
- LPCAT3

Validated extended axes:
- ACSL3/SCD MUFA resistance
- MBOAT1/2 MUFA phospholipid remodeling
- peroxisome/ether-lipid axis (FAR1, AGPS, GNPAT, AGPAT3, TMEM164)
- FADS2/ELOVL5 endogenous PUFA synthesis

### F3 — Lipid-peroxidation machinery
Core question:
> Is there molecular capacity to drive lethal phospholipid oxidation?

Core:
- POR
- CYB5R1

Extended:
- ALOX15
- PHKG2

Protective detoxification:
- PLA2G6

PEBP1 is mechanistic because its relevant function is complex/scaffold dependent.

### F4 — Anti-ferroptotic defense
Core question:
> Are the major lipid-antioxidant systems sufficiently active to suppress ferroptosis?

Core axes:
- System Xc− / GPX4
- FSP1 / CoQ
- DHODH / mitochondrial CoQ
- GCH1 / BH4

Extended:
- GSH synthesis machinery
- transsulfuration
- DHFR/BH4 regeneration

### F5 — Context / response
Examples:
- NFE2L2/KEAP1
- HMOX1
- CHAC1
- SAT1

These genes can be biologically informative but must not dominate the mechanistic score.

## RNA interpretability rule

A gene can be mechanistically important and still be excluded from direct RNA scoring.

Examples:
- NFE2L2 activity depends heavily on KEAP1-mediated protein stabilization.
- KEAP1 redox sensing is not captured by transcript abundance.
- ATG5/ATG7 are broad autophagy genes.
- PEBP1 acts through a protein complex with lipoxygenase.
- HMOX1 can have opposing ferroptosis effects depending on context.

## Core-vs-extended principle

The initial core is intentionally small.

Current provisional core candidates:
- NCOA4
- TFRC
- ACSL4
- LPCAT3
- POR
- CYB5R1
- SLC7A11
- GPX4
- AIFM2/FSP1
- DHODH
- GCH1

Everything else must earn additional weight through empirical reproducibility.

## Important caveat

Even the core is **not yet frozen**.

The next step is to test each RNA-eligible gene across:
- rescue-validated RSL3 datasets;
- cross-inducer Erastin/RSL3 datasets;
- held-out orthogonal induction;
- generic oxidative stress controls;
- neuronal ferroptosis datasets.

A mechanistically strong gene may receive low final RNA weight if its transcript is not reproducible or informative.
