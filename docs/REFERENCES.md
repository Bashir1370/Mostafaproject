# Primary References and Reference Datasets

## Necroptosis mechanistic foundation

1. **Sun L, Wang H, Wang Z, et al.** Mixed lineage kinase domain-like protein mediates necrosis signaling downstream of RIP3 kinase. *Cell*. 2012;148:213–227.  
   DOI: https://doi.org/10.1016/j.cell.2011.11.031  
   Role: identifies MLKL as a functional RIPK3 substrate and essential downstream necroptosis mediator.

2. **Cai Z, Jitkaew S, Zhao J, et al.** Plasma membrane translocation of trimerized MLKL protein is required for TNF-induced necroptosis. *Nature Cell Biology*. 2014.  
   Role: MLKL oligomerization/translocation as terminal execution biology.

3. **Hildebrand JM, Tanzer MC, Lucet IS, et al.** Activation of the pseudokinase MLKL unleashes the four-helix bundle domain to induce membrane localization and necroptotic cell death. *PNAS*. 2014;111:15072–15077.  
   DOI: https://doi.org/10.1073/pnas.1408987111

## Necroptosis relevance to chemotherapy-induced neuropathic pain

4. **RIP3/MLKL pathway-regulated necroptosis: A new mechanism of paclitaxel-induced peripheral neuropathy.** *Journal of Biochemical and Molecular Toxicology*. 2021.  
   PMID: 34056794  
   DOI: https://doi.org/10.1002/jbt.22834  
   Role: direct experimental evidence linking RIP3/MLKL-dependent necroptosis to paclitaxel-induced peripheral neuropathy.

5. **Ma D, Wang X, Liu X, et al.** Macrophage Infiltration Initiates RIP3/MLKL-Dependent Necroptosis in Paclitaxel-Induced Neuropathic Pain. *Mediators of Inflammation*. 2022;2022:1567210.  
   PMID: 36164389  
   DOI: https://doi.org/10.1155/2022/1567210  
   Role: TNF-rich macrophage context, pMLKL, DRG neuronal necroptosis and pain behavior.

## Ferroptosis mechanistic foundation

6. **Doll S, Proneth B, Tyurina YY, et al.** ACSL4 dictates ferroptosis sensitivity by shaping cellular lipid composition. *Nature Chemical Biology*. 2017;13:91–98.  
   DOI: https://doi.org/10.1038/nchembio.2239

7. **Gao M, Monian P, Pan Q, et al.** Ferroptosis is an autophagic cell death process. *Cell Research*. 2016;26:1021–1032.  
   DOI: https://doi.org/10.1038/cr.2016.95

8. **Zou Y, Li H, Graham ET, et al.** Cytochrome P450 oxidoreductase contributes to phospholipid peroxidation in ferroptosis. *Nature Chemical Biology*. 2020;16:302–309.  
   DOI: https://doi.org/10.1038/s41589-020-0472-6

9. **Bersuker K, Hendricks JM, Li Z, et al.** The CoQ oxidoreductase FSP1 acts parallel to GPX4 to inhibit ferroptosis. *Nature*. 2019;575:688–692.  
   DOI: https://doi.org/10.1038/s41586-019-1705-2

10. **Doll S, Freitas FP, Shah R, et al.** FSP1 is a glutathione-independent ferroptosis suppressor. *Nature*. 2019;575:693–698.  
    DOI: https://doi.org/10.1038/s41586-019-1707-0

11. **Mao C, Liu X, Zhang Y, et al.** DHODH-mediated ferroptosis defence is a targetable vulnerability in cancer. *Nature*. 2021;593:586–590.  
    DOI: https://doi.org/10.1038/s41586-021-03539-7

12. **Soula M, Weber RA, Zilka O, et al.** Metabolic determinants of cancer cell sensitivity to canonical ferroptosis inducers. *Nature Chemical Biology*. 2020;16:1351–1360.  
    DOI: https://doi.org/10.1038/s41589-020-0613-y

## NEC reference datasets

- **GSE108621** — human HT-29; DMSO, TNF, TSZ, TSZ+Nec-1s; n=3/group. Primary rescue/inflammation-controlled calibration.
- **GSE172027** — human astrocytes; DMSO, Nec-1s, TSZ, TSZ+Nec-1s; n=3/group. Human neural-lineage rescue replication.
- **GSE154230** — mouse microglia/astrocytes; RIPK1-activating stimuli ± Nec-1s. Neuroimmune robustness.
- **GSE134234** — engineered mouse RIPK3 direct dimerization time course with TNF control. Orthogonal validation.
- **GSE268650** — human optogenetic RIPK3 oligomerization with inhibitor/light controls. Orthogonal human validation.
- **GSE287439** — human iPSC neuron/astrocyte/microglia tri-culture; TSZ ± RIPK1 inhibitor; scRNA-seq. Cell-type validation.

## FPT reference datasets

- **GSE182638** — human MM1R/MM1S; untreated, RSL3, RSL3+Fer-1; n=3/condition/cell line; raw counts. Primary rescue calibration.
- **GSE247883** — human A549; DMSO, RSL3, RSL3+Fer-1; n=3/group. Independent rescue replication.
- **GSE255459** — human BT549/HS578/SUM159; DMSO, Erastin, RSL3; duplicates. Cross-inducer consensus.
- **GSE319384** — human PANC-1/TP53-KO; DMSO, Erastin, RSL3, Ferroptocide; n=3. Held-out multi-inducer validation.
- **GSE287284** — mouse neuronal N2a; DMSO vs RSL3; n=3. Neural validation.
- **GSE152988** — human iPSC-derived neuronal genetic ferroptosis model. Human-neural robustness.

## Generic-stress comparator

- **GSE104664** — human HUVEC; untreated vs H2O2 16h/36h; n=3/group; raw counts.

## Reference policy

- Prefer primary mechanistic experiments.
- Prefer induction plus pathway-specific rescue.
- Prefer raw counts or reproducibly reprocessed SRA.
- Distinguish calibration, replication and held-out validation datasets.
- Do not count multiple arms from one publication as independent studies.
- Record every reference-library change in the changelog.
