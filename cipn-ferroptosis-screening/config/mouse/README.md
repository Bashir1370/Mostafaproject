# Locked mouse gene-set layer

This directory is shared by the two mouse DRG validation datasets (`GSE125002`, `GSE286387`) and is deliberately separate from the rat locks used by `GSE160543`.

`vinik_2024_24_human_to_mouse_locked.csv` is the manually reviewed Human->Mouse Vinik mapping. `00_prepare_mouse_gene_sets.R` audits that mapping with babelgene and freezes the four Human-MSigDB->Mouse sets plus their version/provenance manifest. Main analysis scripts refuse to run until those MSigDB lock files exist.

This preserves the same conceptual source definitions across rat and mouse while making the computational ortholog step explicit. These are cross-species mappings, not evidence that the original human signatures were experimentally validated in mouse DRG.

## Audit exceptions and deduplication

The Vinik-24 audit preserves the raw live babelgene result rather than forcing agreement with the locked mapping. For human `GARS1`, babelgene currently returns the legacy mouse symbol `Gars`, whereas NCBI/MGI list `Gars1` as the current official mouse symbol (Gene ID 353172) and `Gars` as an alias. The locked mapping therefore remains `Gars1` and is explicitly labeled `MANUAL_REVIEW_ACCEPTED_SYMBOL_NOMENCLATURE_EXCEPTION`.

For MSigDB ortholog mappings, multiple human members can collapse onto the same mouse symbol. Frozen membership is therefore unique by `signature + mouse_gene_symbol`; when duplicate mapped rows occur, the row with the larger `num_ortholog_sources` is retained. This prevents duplicated genes from being represented more than once in downstream enrichment sets.

