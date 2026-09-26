# Scoring Specification

## Symbols

- `M_g`: mechanistic prior
- `E_g`: empirical reproducibility
- `S_g`: specificity
- `W_MPS,g`: permissiveness weight
- `W_DSI,g`: discrimination weight
- `r_gs`: within-sample expression rank
- `q_gs`: direction-adjusted gene score

## Gene weights

```text
W_MPS,g = M_g × E_g
W_DSI,g = M_g × E_g × S_g
```

## Direction-aware gene score

Promoter:

```text
q_gs = r_gs
```

Protective gene:

```text
q_gs = 1 - r_gs
```

## Stage score

```text
Stage_k(s) =
  Σ(W_g × q_gs) / Σ(W_g)
```

## DPT

```text
DPT_init = sqrt(D1 × D2)
DPT-MPS  = (DPT_init + D3) / 2
DPT-MCI  = (D1 × D2 × D3)^(1/3)
```

- D1 = cystine loading
- D2 = reducing-capacity vulnerability
- D3 = WRC/actin execution permissiveness

## FPT

```text
FPT_init = (F1 × F2 × F3)^(1/3)
FPT-MPS  = (FPT_init + F4) / 2
FPT-MCI  = (F1 × F2 × F3 × F4)^(1/4)
```

- F1 = iron availability
- F2 = PUFA-phospholipid susceptibility
- F3 = lipid-peroxidation machinery
- F4 = anti-ferroptotic defense failure

## Per-study empirical evidence

```text
Magnitude_gj  = percentile_rank(|beta_gj|)
Confidence_gj = 1 - lfsr_gj

e_gj =
  sign(d_g × beta_gj)
  × Magnitude_gj
  × Confidence_gj
```

## Across-study empirical score

```text
Positive_g      = mean(max(e_gj, 0))
Contradiction_g = mean(max(-e_gj, 0))

E_g = max(0, Positive_g - Contradiction_g)
```

## Specificity

```text
C_g = max(
  competing_program_same_direction,
  generic_stress_same_direction
)

S_g = max(
  0,
  (E_target,g - C_g) /
  (E_target,g + C_g + epsilon)
)
```

## DSI

```text
DSI =
  (D_specific - F_specific) /
  (D_specific + F_specific + epsilon)
```

Range: [-1, +1].

**DSI is relative specificity, not probability.**

## Mandatory safeguards

1. Directly manipulated genes cannot validate themselves in the same experiment.
2. Contradictory effects receive negative evidence.
3. Multiple contrasts from one study are aggregated before cross-study combination.
4. Independent studies are the replication unit.
5. Shared redox biology can contribute to permissiveness but must not dominate DSI.
6. Protein-level damage substrates are not assumed to be RNA activity markers.
7. Thresholds are calibrated only on reference data.
