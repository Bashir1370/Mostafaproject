# Mouse oxaliplatin validation datasets

This project now treats each GEO accession as an independent analysis unit. Results must never be mixed across accessions before a formal cross-dataset synthesis.

## GSE125002 -- mouse DRG, short/local oxaliplatin model
- Organism: *Mus musculus* (C57BL/6)
- Tissue: whole L3-L5 dorsal root ganglia
- Contrast used here: Oxaliplatin vs Vehicle only (cisplatin samples excluded)
- RNA-seq: Illumina NextSeq 500, 75-nt single-end; STAR + HTSeq counts deposited in GEO
- Design: n=3 biological replicates/group; each replicate pools DRG from four animals
- Oxaliplatin: 40 ugg cumulative intraplantar
- Purpose here: independent mouse DRG replication of the ferroptosis-associated transcriptomic screen.
- Limitation: small n and pooled animals; treatment route differs from the rat GSE160543 model.

## GSE286387 -- mouse DRG, chronic systemic oxaliplatin model
- Organism: *Mus musculus* (C57BL/6 male)
- Tissue: bilateral lumbar DRG (L1-L6 reported in the source study)
- Contrast: 10 mg/kg oxaliplatin vs vehicle
- Treatment: intraperitoneal, once weekly for 8 weeks
- RNA-seq analysis used in the source study: n=5/group
- Purpose here: chronic, systemic and newer independent mouse DRG validation.
- Source study reports 13,145 detected genes, 820 DEGs at FDR<0.05, 271 at FDR<0.05 and |FC|>=1.5, and validation of TXNIP at RNA/protein level.
- Loader policy: analysis proceeds only if a genome-scale raw-count matrix is recovered. Normalized-only or DEG-only tables are rejected for DESeq2/GSEA.

## Shared gene-set policy
Both mouse datasets use the same five conceptual sets as GSE160543: Vinik-24 (primary), GOBP_FERROPTOSIS and WP_FERROPTOSIS (supporting), HALLMARK_APOPTOSIS and HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY (context controls). To preserve cross-species comparability, the MSigDB source remains the human database and is explicitly ortholog-mapped to mouse, then frozen. Vinik-24 is likewise mapped human->mouse and frozen. Cross-species mapping is computational and is not equivalent to mouse validation of the original human signature.
