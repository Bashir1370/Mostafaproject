# Scoring Specification — NEC–FPT

## Symbols

- `M_g`: mechanistic prior
- `E_g`: empirical reproducibility
- `S_g`: specificity
- `W_MPS,g`: permissiveness weight
- `W_NFSI,g`: discrimination weight
- `r_gs`: within-sample expression rank
- `q_gs`: direction-adjusted gene score

## Gene weights

```text
W_MPS,g  = M_g * E_g
W_NFSI,g = M_g * E_g * S_g
```

## Direction-aware score

Promoter:
```text
q_gs = r_gs
```

Protective:
```text
q_gs = 1-r_gs
```

## Stage score

```text
Stage_k(s) =
  sum(W_g*q_gs) / sum(W_g)
```

## Necroptosis

```text
N1 = RIPK1/RIPK3 necrosome competence
N2 = MLKL execution competence

NEC-MPS = (N1 + N2)/2
NEC-MCI = sqrt(N1*N2)
```

Trigger-context and checkpoint genes are annotations unless empirical calibration supports transcriptomic use.

## Ferroptosis

```text
FPT_init = (F1*F2*F3)^(1/3)
FPT-MPS  = (FPT_init+F4)/2
FPT-MCI  = (F1*F2*F3*F4)^(1/4)
```

- F1 iron availability
- F2 PUFA-phospholipid susceptibility
- F3 lipid-peroxidation machinery
- F4 defense

## Per-contrast empirical evidence

```text
Magnitude_gj  = percentile_rank(|beta_gj|)
Confidence_gj = 1-lfsr_gj

e_gj =
  sign(d_g*beta_gj)
  * Magnitude_gj
  * Confidence_gj
```

## Across-study E

```text
Positive_g      = mean(max(e_study,0))
Contradiction_g = mean(max(-e_study,0))

E_g = max(0, Positive_g-Contradiction_g)
```

## Specificity

```text
C_g = max(
  competing_program_same_direction,
  generic_stress_same_direction
)

S_g = max(
  0,
  (E_target,g-C_g) /
  (E_target,g+C_g+epsilon)
)
```

## NFSI

```text
NFSI =
  (NEC_specific-FPT_specific) /
  (NEC_specific+FPT_specific+epsilon)
```

Range [-1,+1].

## Mandatory safeguards

1. Directly manipulated genes cannot self-validate.
2. Post-translational execution is not inferred from RNA alone.
3. Contradictory effects receive negative evidence.
4. Multiple contrasts/cell lines are aggregated within study.
5. Independent study is the replication unit.
6. Rescue designs are prioritized.
7. TNF-only inflammatory behavior is explicitly controlled in NEC calibration where available.
8. Generic oxidative stress is explicitly controlled for FPT specificity.
9. Thresholds are calibrated only from reference data.
