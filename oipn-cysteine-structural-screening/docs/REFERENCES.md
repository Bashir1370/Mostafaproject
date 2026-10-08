# Reference register

Sources supporting the design; retrieval/annotation/software releases will be locked at implementation. Study references do not imply that all supplementary data have already been downloaded or checked.

| ID | Source | Role / boundary |
|---|---|---|
| R01 | [A structural bioinformatics framework for prioritizing pH-sensitive proteins from 3D structural features](https://doi.org/10.1186/s12859-026-06415-1), 2026 | Motivating workflow; histidine/tumor pH descriptors are not directly transferable |
| R02 | [Prediction of reversibly oxidized protein cysteine thiols using protein structure properties](https://doi.org/10.1110/ps.073252408), Sanchez et al., 2008; PMID 18287280 | Cys accessibility, pKa and nearest-Cys distance; published COPA is not our equal-weight formula |
| R03 | [Functional site profiling and electrostatic analysis of cysteines modifiable to cysteine sulfenic acid](https://doi.org/10.1110/ps.073096508), Salsbury et al., 2008; PMID 18227433 | Local structural context; surface exposure/secondary structure are not universal standalone classifiers |
| R04 | [Systematic and Quantitative Assessment of Hydrogen Peroxide Reactivity With Cysteines Across Human Proteomes](https://doi.org/10.1074/mcp.RA117.000108), Fu et al., 2017; PMID 28827280 | Candidate independent experimental benchmark; human-cell context, not mouse-DRG proof |
| R05 | [Benchmarking In Silico Tools for Cysteine pKa Prediction](https://doi.org/10.1021/acs.jcim.3c00004), 2023; PMID 36996330 | pKa uncertainty and predictor limitations |
| R06 | [AlphaFold DB FAQ](https://www.alphafold.ebi.ac.uk/faq) | Local confidence semantics; project thresholds are pragmatic QC |
| R07 | [Global Inhibition of Reactive Oxygen Species (ROS) Inhibits Paclitaxel-Induced Painful Peripheral Neuropathy](https://doi.org/10.1371/journal.pone.0025212), 2011 | Preclinical support for ROS involvement; different chemotherapy |
| R08 | [Oxidative stress in the development, maintenance and resolution of paclitaxel-induced painful neuropathy](https://doi.org/10.1016/j.neuroscience.2016.06.050), 2016 | ROS findings depend on time, cell type and assay/context |
| DATA01 | [GSE286387](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286387); [source study: Chronic Oxaliplatin Treatment Induces CIPN in Mice via Activation of the TXNIP Pathway](https://pmc.ncbi.nlm.nih.gov/articles/PMC13356913/) | Candidate chronic systemic mouse DRG discovery dataset; current audit pending |
| DATA02 | [GSE125002](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125002) | Candidate complementary short/local mouse DRG dataset; distinct biological design |

Repository context reviewed before registration: root NEC–FPT workflow; ferroptosis-only workspace; cipn-ferroptosis-screening mouse metadata, loaders and shared functions. They are not validation datasets for this structural index merely because they contain ROS/ferroptosis signatures.
