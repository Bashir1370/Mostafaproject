# Changelog

## v0.2.0 — Necroptosis–Ferroptosis redesign

### Breaking changes
- Replaced **disulfidptosis (DPT)** with **necroptosis (NEC)** as the second active regulated-cell-death program.
- Replaced `DSI` with **NFSI (Necroptosis–Ferroptosis Specificity Index)**.
- Removed the DPT gene dictionary from the active repository.
- Rebuilt reference datasets and contrasts around NEC/FPT perturbation studies.

### Added
- NEC core architecture:
  - N1 RIPK1/RIPK3 necrosome competence
  - N2 MLKL execution competence
- NEC trigger/checkpoint annotation layer.
- GSE108621 primary NEC rescue/inflammation-controlled calibration.
- GSE172027 human astrocyte rescue replication.
- GSE154230 neuroimmune robustness.
- GSE134234 and GSE268650 orthogonal RIPK3 validation.
- GSE287439 human neural tri-culture scRNA validation.
- GSE182638 primary FPT rescue calibration.
- GSE319384 held-out multi-inducer FPT validation.
- neural FPT validation datasets.

### Preserved
- mechanistic prior M
- empirical reproducibility E
- specificity S
- direction-aware evidence
- anti-circularity
- study-level replication rule
- MPS/MCI/ESR separation
- bootstrap / whole-study validation policy

### Next
- run the updated reference inspection;
- calibrate GSE108621 and GSE182638 first;
- then expand to replication datasets.

## v0.1.0 — Historical DPT–FPT architecture
The original DPT–FPT framework remains available through Git history but is no longer the active design.
