# Step 05b2 — hash-gated canonical cysteine mapping

Protocol 0.1.0 is unchanged. This mapper finds canonical Cys/residue/SG correspondences in every metadata-eligible deposited chain/model. It does not select a best conformation or approve local quality, chemical state, biological assembly or oxidation susceptibility.

## Dependencies and Ubuntu execution

Use an isolated subproject environment and the pinned Gemmi 0.7.5 binary wheel for genuine mmCIF parsing. The release provides CPython 3.13 / Linux x86_64 manylinux wheels compatible with the workstation's glibc; the validation runtime used its Python 3.11 wheel. No R/global Python packages are changed.

```bash
cd /home/bashir/Desktop/Mostafaproject
git pull --ff-only
python3 -m venv oipn-cysteine-structural-screening/.venv
oipn-cysteine-structural-screening/.venv/bin/python -m pip install --only-binary=:all: -r oipn-cysteine-structural-screening/requirements-structure.txt
oipn-cysteine-structural-screening/.venv/bin/python oipn-cysteine-structural-screening/scripts/05b2_map_cysteines_GSE286387.py
cat oipn-cysteine-structural-screening/results/05b2_cysteine_mapping/mapping_report.md
```

After dependency installation, mapping is entirely offline and uses existing raw caches. Keep original data/raw/structure_coordinates/step05b1/ files and sidecars. No structural downloads are repeated. Inputs are frozen in config/GSE286387_step05b2_input.json; every bound upstream input, recorded acquisition input and all 2,638 raw file/sidecar pairs must pass accepted byte-hash/size/provenance checks before parsing begins. A missing/changed raw file fails execution rather than silently becoming an absent residue. Parser/schema execution failures invalidate SUCCESS.txt and produce FAILURE.txt. Expected unresolved scientific correspondence is explicitly held, not an execution exception.

## Correspondence rules fixed before workstation inspection

1. Resolve PDBe candidate author-chain IDs to label subchains using explicit scheme/atom identifiers. Ignore nonpolymer ligand/water subchains sharing an author ID; retain polymer subchains and every deposited model. Keep label_seq_id, auth_seq_id and insertion codes separate. Missing residues are inferred from deposited polymer/reference correspondence, not an observed-atom-only sequence.
2. For AlphaFold, require contiguous deposited polymer sequence IDs and exact sequence identity to the frozen canonical interval, including partial-model boundaries.
3. For PDB, prefer exact deposited UNP accession/entity/chain reference segments. Segment endpoints must lie in the reference, have equal lengths, and have no unresolved insertion codes. Every residue must match canonical sequence, except a substitution explicitly documented at both the exact deposited and database positions with matching monomer identities. Overlapping contradictory/nonunique segments are held. Gapped segments, unexplained mismatches or unresolved deposited UNP associations are held for later explicit per-residue evidence review; no heuristic alignment or rescue by homolog name is used.
4. Only when no deposited UNP entity association exists may a complete, contiguous, known deposited polymer sequence map by one unique exact occurrence in the frozen canonical sequence. Repeated/missing matches or unknown residues are held. No similarity cutoff or score-driven choice is introduced.
5. For each canonical Cys, distinguish outside-checked-segment, residue unobserved, substituted/modified target, missing SG, unresolved author identifiers and mapped native CYS SG. Structural modifications/mutations remain evidence; no atom or chemical state is repaired. Microheterogeneous polymer positions are held.
6. SG must have declared sulfur element and finite coordinates. Select highest known positive occupancy, then blank-before-lexical alternate ID; duplicate SG alternate identifiers or unknown occupancies are held. This chooses an SG record only: other atoms/conformers are not prepared. missing_CYS_heavy_atoms reports absent N/CA/C/O/CB/SG names among deposited target atom records; it is not proof of a complete coherent alternate conformer. SG B_iso_or_equiv is recorded raw and is not used as a positive score or quality pass.

Neither the deposited reference nor a matched SG alone establishes native context. Declared changed canonical positions are exported; neighborhood mutation checks, assembly partners, boundary/context uncertainty and experimental local validation remain subsequent eligibility work. No pLDDT thresholds are applied in this substep; the frozen target >=90 / modeled-neighbor >=70 within 6 A rule remains future local-quality work.

## Outputs and complete accounting

- site_structure_mapping.csv: one record per canonical Cys/candidate/resolved label-subchain/deposited model (held chain correspondences use explicit unresolved model fields). It preserves identities, numbering, component/observation flags, SG selection/coordinates/raw B value, atom-name gaps, declared sequence changes and exact raw-file provenance.
- candidate_mapping_audit.csv: every one of the 1,872 offered catalogue records, including 324 metadata held and 15 excluded records that are not parsed as eligible structures.
- site_step_audit.csv: every one of the 10,799 unique canonical Cys sites, with mapping-option counts and held/unassessable status. Mapped sites remain held for quality review; metadata ambiguity is retained even without an eligible structure.
- protein_mapping_coverage.csv: all 760 proteins, including the 29 without metadata-eligible candidates; coverage is descriptive.
- raw_file_checksums.csv, input_checksums.csv, summary.json, mapping_report.md and generation markers.

Every structural_eligibility/local_quality_status remains not_assessed. SUCCESS.txt means STEP05B2_GENERATED_REVIEW_PENDING, not final Step 05 acceptance. Review held reasons and actual coverage before proceeding. No zero susceptibility, SASA, pKa or score is assigned.

## Validation and practical limits

The 73-test full regression suite passed; all 14 mapping tests, including an additional end-to-end sample integration, passed. Cases include insertion/author versus label numbering, unique and repeated sequence placement, explicit reference offsets, gapped/unexplained references, documented substitutions, missing residues/SG, all-model preservation, altloc occupancy/ties, zero occupancy, microheterogeneity, nonpolymer-chain sharing, metadata ambiguity and local raw-hash rejection. The end-to-end representative fixture retains the complete 760-protein / 10,799-site audit while processing two accepted-hash files.

Accepted-byte real examples: 11gl chains A/B map four SG sites each among eleven canonical Cys, with seven outside the checked fragment; AF-A0A087WRH0-F1 maps all 35 Cys/SG sites. Combined unique mapped sites in the sample are 39. These are validation samples, not the user's complete 1,176-file run. Full corpus execution/runtime/mapping coverage remains pending on Ubuntu; the validator has only representative raw files locally. Strict holds are expected and must not be relaxed to achieve a desired number of structures. Shared per-file atom indexing avoids scanning unrelated polymer/ligand rows for every candidate chain.

## Primary documentation

- Gemmi parser/API: https://gemmi.readthedocs.io/en/stable/cif.html
- Gemmi pinned distribution/wheels: https://pypi.org/project/gemmi/0.7.5/
- seq_align_beg is a pointer to entity_poly_seq.num, not author numbering: https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_struct_ref_seq.seq_align_beg.html
- Label/author/sequence scheme: https://mmcif.pdb.org/docs/user-guide/guide.html
