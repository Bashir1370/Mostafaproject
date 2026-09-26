# Changelog

## v0.1.0 — Architecture freeze candidate

### Added
- Disease-agnostic DPT-vs-FPT classifier scope.
- DPT stages:
  - cystine loading
  - reducing-capacity vulnerability
  - WRC/actin execution
  - empirical state resemblance
- FPT stages:
  - iron availability
  - PUFA-phospholipid susceptibility
  - lipid-peroxidation machinery
  - anti-ferroptotic defense failure
- Mechanistic evidence tiers A+/A/B/C.
- Direction-aware RNA scoring.
- Separate MPS, MCI, ESR and DSI concepts.
- DPT cystine × reducing-capacity interaction gate.
- Mechanistic prior M, empirical reproducibility E and specificity S.
- Anti-circularity rule for directly manipulated genes.
- Contradiction penalty across studies.
- Generic oxidative-stress comparator strategy.
- Four-class output:
  - DPT-dominant
  - FPT-dominant
  - Mixed
  - Indeterminate / Neither
- Reference-data-based threshold calibration.
- Initial DPT and FPT gene dictionaries.
- Initial perturbation dataset registry.
- Primary reference bibliography.

### Not yet completed
- Reference-dataset preprocessing.
- Differential-expression re-analysis.
- Empirical E calculation.
- Specificity S calculation.
- Final MPS/DSI weights.
- Threshold calibration.
- Leave-one-study-out validation.
- Software implementation and unit tests.
