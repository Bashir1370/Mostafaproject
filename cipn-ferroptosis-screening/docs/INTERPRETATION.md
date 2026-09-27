# Interpretation guide

## What this analysis can answer

The analysis can support statements such as:

> Oxaliplatin-treated DRG shows significant enrichment of independent ferroptosis-associated transcriptomic signatures.

or:

> No transcriptomic evidence supporting ferroptosis-associated activation was detected in this DRG dataset.

## What it cannot prove

Bulk DRG RNA-seq cannot, by itself, establish:
- lipid peroxidation,
- GPX4 enzymatic inhibition,
- Fe2+ accumulation,
- ferroptotic membrane damage,
- or the identity of the DRG cell type producing the signal.

Therefore avoid writing simply:

> Ferroptosis occurred.

Prefer:

> The DRG transcriptome showed strong/suggestive/no evidence of a ferroptosis-associated transcriptional state.

## Why the Vinik-24 set is primary

KEGG and WikiPathways contain mechanistic components with mixed biological directions. For example, anti-ferroptotic defense genes may increase as a compensatory response during ferroptotic stress.

The Vinik-24 panel is used as the primary directional transcriptomic signature because it was selected for ferroptosis-vs-apoptosis discrimination and validated as an induced biomarker panel.

## Why apoptosis and ROS are included

CIPN causes broad injury and oxidative stress. A positive ferroptosis pathway alone is not automatically specific.

The control signatures help identify cases where the transcriptome reflects a general injury/death response rather than a selective ferroptosis-like state.

## Suggested manuscript language

### Strong support
"Multiple independent ferroptosis-associated gene sets were positively enriched in treated DRG, including the experimentally validated Vinik-24 transcriptional biomarker panel and a curated ferroptosis pathway, providing strong transcriptomic evidence consistent with ferroptosis-associated activation."

### Suggestive
"Ferroptosis-associated signatures showed positive but incomplete convergence, supporting a possible ferroptosis-related transcriptional response that requires orthogonal validation."

### No support
"No consistent enrichment of independent ferroptosis-associated transcriptomic signatures was observed."

## Experimental confirmation, if available

The computational result becomes much stronger if supported by one or more orthogonal assays:
- C11-BODIPY lipid ROS
- 4-HNE or MDA
- labile Fe2+
- GPX4 protein/activity
- ACSL4 protein
- Ferrostatin-1 or Liproxstatin-1 rescue
