# Ten-step protocol — v0.1.0

## Scope and registration

Species: Mus musculus (mouse), tissue: DRG, treatment: oxaliplatin, control: matched vehicle. Discovery input is significant protein-coding DEGs, BH FDR < 0.05, both directions, no fold-change cutoff. Dataset choice is audited before differential expression. The structural score is fixed before inspection of favored candidates. v0.1.0 freezes the design, not an empirically validated predictor.

There are **exactly ten scientific steps**. Environment preparation is operational setup, not an additional scientific step. Every step emits an audit with object ID, status, reason, counts entering/leaving, input checksums, software versions and protocol version. A rejected site need not reject its protein. Unassessable is never encoded as zero susceptibility.

## 01 — Dataset selection and sample QC

**Input:** GEO records, source publication, raw/processed data and metadata.

**Entry:** mouse DRG; oxaliplatin versus matched vehicle; at least three independent biological RNA-seq units per group; genome-scale expression data; documented treatment and collection time; neuropathy phenotype supported by the source study. A pooled library counts as one biological unit, not as the number of pooled mice.

**Procedure:** verify species, DRG levels, sex, age, strain, route, dose, schedule, collection time, pooling, animal IDs and batch. Match every column to exactly one audited sample. Inspect library size, detected genes, sample correlation and PCA. Raw gene counts or raw reads reprocessed into counts are required for the chosen DESeq2 route; apparent integer values alone do not prove raw-count provenance. Never round TPM/FPKM into counts. No sample is excluded merely because it is distant in PCA; exclusions require documented technical/identity evidence.

**Exit/hold:** wrong species/tissue; inseparable combination therapy; unresolved identity; fewer than three independent units per group; complete treatment–batch confounding; only DEG/normalized tables available and no recoverable count route. Stop inference while a count/sample audit is unresolved. Do not switch datasets based on how many desired genes are found.

**Output:** audited sample manifest, data/source manifest, complete count matrix and QC report. Candidate primary: GSE286387 (chronic systemic model); complementary: GSE125002 (different short/local model). Current acceptance (2026-10-08): GSE286387 is the primary dataset, all ten samples retained, unpaired `~ condition`, oxaliplatin/control. See [design evidence and limitations](DATASET_DESIGN_AUDIT.md). GSE125002 remains unaudited and separate.

## 02 — Differential expression and discovery gene list

**Input:** accepted counts and sample design.

**Entry:** retain a gene if raw count >= 10 in at least min(n_vehicle, n_oxaliplatin) samples across the two groups. Preserve stable mouse gene IDs and annotation release. Do not merge different stable genes merely because symbols coincide.

**Procedure:** DESeq2; numerator Oxaliplatin, denominator Vehicle. Use design ~ condition unless the audited experiment requires estimable batch/block terms; record the final model before fitting. Use independentFiltering=FALSE to apply the declared count filter consistently, while retaining standard DESeq2 outlier handling. BH correction is over the valid p-values of all tested genes, not just selected protein-coding candidates. Keep all results and NA reasons.

**Pass to step 03:** protein-coding gene, finite padj < 0.05. No absolute log2FC or baseMean ranking threshold. Both up- and down-regulation enter.

**Exit:** failed count filter; noncoding biotype; missing/invalid adjusted p-value; padj >= 0.05. Non-significance is not evidence of oxidative resistance. Zero DEGs produces a documented empty discovery result; never relax FDR to obtain candidates.

**Output:** all tested genes, significant protein-coding DEGs, retained tested-gene universe and exclusion audit.

## 03 — Gene-to-protein and sequence mapping

**Input:** significant protein-coding mouse genes.

**Entry:** documented gene-to-protein mapping in Mus musculus and retrievable sequence. Use the UniProt-designated canonical sequence when available. If no canonical mapping resolves ambiguity, hold the gene for review rather than choosing the longest or first arbitrary match. Record database release/access date and sequence checksum.

**Exit/hold:** missing sequence, wrong species, unresolved gene/protein identity, or unresolved canonical mapping. RNA at gene level does not demonstrate which isoform is expressed. Isoform alternatives are secondary analyses.

**Output:** gene–protein map, canonical FASTA, sequence length, species and mapping audit. The executable first pass is specified in [STEP_03.md](STEP_03.md): exact current Ensembl GeneId links, all mouse entry candidates retained, multiple accessions held, displayed-sequence association checked, and ambiguous/fragment sequences unassessable. Historical IDs without current links require separate rescue review. Canonical means the UniProtKB representative sequence, not a measured dominant DRG isoform. The subsequent [Step 03b](STEP_03B.md) resolves multiple entries only with an explicit gene-centric representative, complete candidate accounting and existing exact-gene mouse sequence evidence. Preserve held cases and report gene coverage; the user requires most genes mapped before Step 04, without forced identities.

## 04 — Cysteine inventory

**Input:** accepted sequences.

**Entry:** at least one cysteine (C). Use 1-based canonical sequence numbering.

**Exit:** no cysteine. Ambiguous/invalid sequence sections are reviewed and not silently treated as resolved amino acids.

**Output:** one row per Cys, total Cys count, Cys percentage and protein length. Counts/density are descriptive, not score inputs. The executable implementation and reviewed input binding are documented in [STEP_04.md](STEP_04.md); U is not counted as C, precursor numbering is preserved, and unresolved upstream genes stay in the full-universe audit.

## 05 — Structure retrieval, mapping and local quality

**Input:** canonical sequences and Cys inventory. The reviewed Ubuntu inventory is bound in config/GSE286387_step05_input.json. [Step 05a](STEP_05A.md) implements a candidate-metadata catalogue; [05b1](STEP_05B1.md) implements raw acquisition and [05b2](STEP_05B2.md) implements canonical Cys correspondence. Local-quality/context acceptance remains separate work within Step 05b. Metadata candidates are not structural passes.

**Priority:** suitable experimental mouse structure first; otherwise predicted mouse structure. Human/rat ortholog structures do not silently substitute in the primary analysis. Resolve isoform and residue numbering by explicit alignment. Select structures by mapping, local completeness, native sequence/context and quality before scoring; never choose the highest-scoring conformation.

**Entry:** mapped Cys present with SG and the atoms needed for SASA/pKa; no engineered mutation at the target or within the local neighborhood; locally interpretable coordinates. For alternate positions, use highest occupancy and a deterministic alternate-ID tie rule. Preserve relevant assembly partners, ligands and cofactors and document preparation. A monomer with unknown assembly context is labeled accordingly.

**Predicted-model rule:** target Cys pLDDT >= 90; all modeled neighboring residues with a non-hydrogen atom within 6 Å of target SG have pLDDT >= 70. The target is excluded from the neighbor set. Missing confidence is not a pass. These are pragmatic QC thresholds, not oxidation cutoffs. Check domain/assembly placement uncertainty where it affects the local environment.

**Experimental rule:** inspect deposited local validation and atom completeness; unresolved local validation/context is held for manual review. Global resolution alone does not replace local assessment. The method-specific implementation of this review must be documented before bulk feature extraction.

**Exit:** target absent, SG absent, wrong mapping, low local confidence, unresolved missing neighborhood, target mutated, or unresolvable context. Partial structures remain eligible for covered sites; report unobserved Cys separately. A fragment may have an artificial exposed surface, so review boundaries near the target.

**Output:** structure manifest, mapping, quality mask, selected model/chain and assembly provenance. A protein with no eligible site is structurally unassessable, not assigned zero.

## 06 — Chemical-state classification

**Input:** eligible structural sites, sequence features and curated/deposited annotations.

**Classes:** compatible_with_free_thiol; structural_disulfide; metal_coordinating; incompatible_modified; unresolved_state. Separate regulatory/reversible disulfide annotations from permanent structural bonds.

**Pass:** compatible_with_free_thiol after the state checks, with annotation completeness recorded. Absence of annotation is not experimental proof of a free thiol.

**Separate from primary score:** permanent structural disulfides, metal ligands and modifications incompatible with the free-thiol model. Hold contradictory/ambiguous states. Do not force a disulfide to reduced state merely to obtain a pKa. Known oxidized/regulatory sites can supply supporting evidence or use a documented reduced-state structure; they are not erased from the evidence registry.

**Output:** state table, evidence source and per-site decision. Retain a protein if it has another eligible site.

## 07 — Feature extraction

**Input:** quality/state-eligible sites and documented prepared structures.

**Features:** solvent-accessible surface area of target SG in Å²; predicted Cys pKa; shortest SG–SG distance to another represented cysteine, excluding itself. Use a consistent SASA implementation, atomic-radii convention and 1.4 Å probe. PROPKA is the planned pKa engine; exact version, invocation and preparation are recorded at implementation. No neighboring Cys gives distance=NA with reason no_other_represented_cys, not zero.

**Pass:** finite nonnegative SG SASA and finite successfully assigned pKa for the mapped site. Zero SASA is valid and is not a hard exclusion. Distance is supplemental and can be missing without blocking the score.

**Exit/hold:** missing primary feature, calculation failure, ambiguous residue identity or unsupported chemical state. Never substitute a default pKa or missing SASA with zero.

**Output:** raw site features, calculation logs and feature-completeness audit. Keep pKa as a prediction; it is not a measured kinetic rate or redox potential.

## 08 — Site scoring, protein ranking and shortlist

**Input:** all sites with both primary features; the score reference population is this full discovery site set.

**Rule:** score = 0.5 * percentile_rank(SG_SASA) + 0.5 * percentile_rank(-pKa). Protein score is maximum eligible-site score. Detailed rank/tie/missing rules are in SCORING_SPEC.md.

**Shortlist:** scores at or above the empirical 75th percentile of scoreable proteins, preserving ties. This is a relative top-quartile rule, not a biochemical susceptibility threshold. Export the full ranking as well. If too few proteins exist for a useful quartile comparison, report an exploratory list and do not present enrichment as robust.

**Exit:** no complete eligible site. Tied scores receive tied scientific ranks; a stable identifier can order display only. Report best Cys, total/evaluated Cys, coverage and structure source.

**Output:** site scores, complete protein ranking, shortlist and exact reference membership. Do not include fold change, FDR, pathway membership, known function, helices, conservation or pLDDT as positive score terms.

## 09 — Independent evaluation and robustness

**Input:** frozen rule, independent experimental site data and alternative eligible structures/QC scenarios.

**Procedure:** compare with SASA-only and pKa-only rankings; matched Cys-count/length controls where appropriate. See VALIDATION_PLAN.md. Do not optimize v0.1.0 weights using OIPN candidate identities. Evaluate sensitivity to weights, quality, structure source/conformation and protein aggregation.

**Exit classification:** benchmark incomplete = unvalidated; no useful improvement/unstable = exploratory structural index; useful independent support = supported for that benchmark scope only. Negative validation is retained. No fixed retrospective performance cutoff is invented after viewing results.

**Output:** benchmark provenance/results, robustness table, limitations and explicit interpretation status. This is an evaluation gate, not proof of in-vivo OIPN oxidation.

## 10 — Biological interpretation

**Input:** full ranking, shortlist, robustness and transcriptomic results.

**Procedure:** annotate functions, location, potential cell types, known redox sites, log2FC/FDR and independent expression evidence. Whole-DRG RNA does not assign a response to neurons rather than glia/immune cells. New/poorly annotated candidates remain reportable. Assess enrichment of the shortlist against all scoreable discovery DEG proteins; use stable gene IDs, Fisher/hypergeometric tests and BH FDR < 0.05 over a predeclared term collection.

**Exit:** no removal merely because a protein is unfamiliar or fails to be DEG in a different model. Lack of replication is recorded, with route/time differences. Do not merge GSE125002 and GSE286387 count matrices or select whichever produces a preferred result.

**Output:** final candidate table, enrichment, selected site explanations, cohort-specific interpretation and limitations. Any claim of actual oxidation, damage or causal pain requires additional direct evidence.

## Audit and amendment rules

For each step record n_input, n_pass, n_excluded, n_held and n_unassessable, with disjoint statuses at the relevant object level. Site and protein counts are separate. Freeze selected dataset, sample list, reference releases, tool versions and structure-selection implementation before the corresponding inference. Amendments receive a version and a reason; v0.1.0 outputs are not overwritten by a redesigned score.
