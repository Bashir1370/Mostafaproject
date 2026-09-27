# Critical Review of Provisional Ferroptosis CORE Genes

## Why this review was necessary

A ferroptosis transcriptomic model should not assume that every well-known ferroptosis regulator is:
1. universally required across all ferroptosis inducers;
2. directionally interpretable at the mRNA level;
3. suitable for the same weight in a sample-level score.

A major cross-screen analysis showed that ferroptosis genetics are highly context-dependent and that different inducers converge more reliably on **lethal plasma-membrane phospholipid peroxidation** than on one universal gene dependency.

Therefore, this project separates:

- **mechanistic strength**
- **universality across induction contexts**
- **RNA interpretability**
- **final empirical transcriptomic usefulness**

The last item will be learned from GEO perturbation datasets rather than assumed from literature.

## Review decisions

### NCOA4 — EXTENDED_ANCHOR
**Keep in F1, but not as a universal CORE gene.**

Why:
- strong causal ferritinophagy/iron evidence;
- loss of NCOA4 can suppress ferroptosis in relevant models;
- however, ferroptosis can obtain redox-active iron through routes other than NCOA4-dependent ferritinophagy.

Interpretation:
NCOA4 is a strong **iron-module anchor**, not a mandatory transcriptomic ferroptosis marker.

### TFRC — EXTENDED
**Demote from CORE.**

Why:
- transferrin receptor biology supports iron availability;
- protein-level TFRC has been reported as a useful ferroptosis marker;
- TFRC transcription is also strongly influenced by proliferation, iron demand and broader iron homeostasis.

Interpretation:
retain for F1, but empirical calibration must determine its RNA weight.

### ACSL4 — CORE
**Retain.**

Why:
- multiple genetic studies show strong control of PUFA-phospholipid susceptibility;
- expression can predict ferroptosis sensitivity in several models;
- however, ACSL4 dependence is not universal and is stronger for some GPX4-inhibition contexts than for cystine deprivation.

Interpretation:
CORE for the lipid-susceptibility module, but never interpreted as proof of ferroptosis by itself.

### LPCAT3 — CORE
**Retain provisionally.**

Why:
- direct genetic evidence supports incorporation/remodeling of oxidizable PUFA phospholipids;
- transcript abundance is more straightforward to interpret than post-translational activity-only nodes.

Interpretation:
retain in the core, subject to cross-study empirical reproducibility.

### POR — CORE
**Retain as a peroxidation anchor.**

Why:
- CRISPR/genetic and lipidomic data directly connect POR to ferroptotic phospholipid oxidation;
- expression is biologically interpretable as oxidative capacity.

Caveat:
POR is not assumed to be the only route to phospholipid peroxidation.

### CYB5R1 — EXTENDED
**Demote from CORE despite strong direct evidence.**

Why:
- knockout evidence directly supports lipid-peroxidation execution;
- the evidence base is narrower across independent systems than ACSL4/GPX4/POR.

Interpretation:
high-priority EXTENDED gene that can return to the final core if GEO reproducibility is strong.

### SLC7A11 — EXTENDED
**Demote from universal CORE.**

Why:
- central suppressor of ferroptosis induced by cystine deprivation/System Xc− inhibition;
- direct GPX4 inhibitors such as RSL3 bypass the upstream requirement for System Xc−;
- stress-responsive transcription can also make directionality less simple in an actively dying sample.

Interpretation:
important defense/permissiveness gene, especially for class-I inducer biology, but not a universal ferroptosis-core transcript.

### GPX4 — CORE
**Retain as the strongest mechanistic defense anchor.**

Why:
- GPX4 is a canonical phospholipid hydroperoxide detoxification system;
- direct GPX4 inhibition or genetic loss can induce ferroptosis.

Critical RNA caveat:
RSL3 inhibits GPX4 activity without requiring GPX4 mRNA to decrease.

Therefore:
GPX4 transcript reflects **defense capacity/permissiveness**, not direct proof of pathway execution.

### AIFM2 / FSP1 — EXTENDED
**Demote from universal CORE but keep high priority.**

Why:
- two independent studies established a GPX4-independent FSP1-CoQ defense axis;
- its quantitative importance varies substantially between cell states.

Interpretation:
strong independent defense module; empirical expression performance determines final weight.

### DHODH — EXTENDED
**Demote from universal CORE.**

Why:
- strong causal mitochondrial ferroptosis-defense biology;
- the effect is particularly pronounced in mitochondrial / GPX4-low contexts.

Interpretation:
high-confidence context-specific defense rather than a universal whole-cell ferroptosis marker.

### GCH1 — EXTENDED
**Demote from universal CORE.**

Why:
- genome-wide activation screening plus metabolic/lipid evidence establish the GCH1-BH4 axis;
- expression can stratify sensitivity in some systems;
- it is nevertheless one of several parallel defense systems rather than a universal requirement.

## Revised provisional CORE

After literature-only review, the **provisional transcriptomic CORE** is intentionally small:

- ACSL4
- LPCAT3
- POR
- GPX4

This does **not** mean these four genes define ferroptosis in isolation.

It means they currently provide the strongest combination of:
- direct mechanistic evidence;
- pathway centrality;
- interpretable transcript abundance;
- and potential cross-context usefulness.

## High-priority EXTENDED anchors

The following must remain in empirical testing and can be promoted after GEO calibration:

- NCOA4
- TFRC
- CYB5R1
- SLC7A11
- AIFM2/FSP1
- DHODH
- GCH1

## Important conceptual conclusion

The model should ultimately be **module-centric, not single-gene-centric**.

Ferroptosis is better represented as the conjunction of:

1. iron availability;
2. oxidizable phospholipid substrate availability;
3. phospholipid-peroxidation capacity;
4. insufficient anti-ferroptotic defense;

rather than as a requirement that one fixed list of transcripts all move in a predetermined direction.

Final RNA weights remain unfrozen until rescue-validated and cross-inducer GEO analyses are completed.
