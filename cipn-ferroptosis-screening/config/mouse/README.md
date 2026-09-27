# Locked mouse gene-set layer

This directory is shared by the two mouse DRG validation datasets (`GSE125002`, `GSE286387`) and is deliberately separate from the rat locks used by `GSE160543`.

`vinik_2024_24_human_to_mouse_locked.csv` is the manually reviewed Human->Mouse Vinik mapping. `00_prepare_mouse_gene_sets.R` audits that mapping with babelgene and freezes the four Human-MSigDB->Mouse sets plus their version/provenance manifest. Main analysis scripts refuse to run until those MSigDB lock files exist.

This preserves the same conceptual source definitions across rat and mouse while making the computational ortholog step explicit. These are cross-species mappings, not evidence that the original human signatures were experimentally validated in mouse DRG.
