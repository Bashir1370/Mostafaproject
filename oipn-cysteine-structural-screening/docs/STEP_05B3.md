# Step 05b3 — local predicted confidence and coherent target atoms

This substep implements the already frozen predicted local-confidence gate. It is an intermediate Step 05 output, not final structural inclusion. Protocol 0.1.0, scientific thresholds, scoring weights and existing executables are unchanged.

## Input and Ubuntu execution

The actual Ubuntu 05b2 mapping is reviewed at source commit 7dd9d5dde4e2a887a3af46f5c6e0c7b90889641c and frozen in docs/audits/ubuntu_cysteine_mapping_2026-10-09. It retains 760 proteins / 10,799 sites / 25,026 mapping records, including 14,310 mapped options. config/GSE286387_step05b3_input.json binds all snapshot hashes and the preceding raw-acquisition/mapping input chain. Every accepted raw body/sidecar is rechecked locally before parsing. Raw cache files must stay at the existing project-relative paths. No network request or new dependency is needed; use the existing pinned Gemmi 0.7.5 environment.

```bash
cd /home/bashir/Desktop/Mostafaproject
git pull --ff-only
oipn-cysteine-structural-screening/.venv/bin/python oipn-cysteine-structural-screening/scripts/05b3_local_quality_GSE286387.py
cat oipn-cysteine-structural-screening/results/05b3_local_quality/quality_report.md
```

## Fixed entry and local gate

Only mapped native CYS/SG options enter this per-option assessment. The entire canonical universe is retained separately, including 265 mapping-unresolved and 793 no-metadata-eligible-structure sites. All 4,569 mapped experimental options are exported as held for experimental local validation/context. Their B factors are never interpreted as pLDDT; global resolution alone does not approve them.

For AlphaFold options:

1. Recheck deposited polymer sequence against the frozen canonical interval. Require one polymer subchain/model, known coordinates, no alternate atoms, unit occupancy, unique atom names and an observed record for every modeled polymer residue. Nonpolymer/ambiguous model context is held; missing atoms are not generated.
2. Require confidence JSON residueNumber to equal model-relative 1..N and finite numeric pLDDT values in [0,100]. Check every CIF atom B_iso_or_equiv against its JSON residue pLDDT within 0.011, accommodating two-decimal serialization. This is a consistency tolerance, not a biological threshold. Unknown/inconsistent source confidence is held.
3. Require the target CYS N/CA/C/O/CB/SG names in the same deposited conformer. In this predicted-model path all alternate IDs must be blank and occupancy one. Confirm the selected sulfur ID/coordinates and canonical position against the accepted mapping; an inconsistent mapping/raw association fails execution.
4. Enumerate all modeled non-H/D atoms within inclusive distance <=6.0 A of target SG in the model. Deduplicate their residue IDs and exclude all atoms of the target residue. Target pLDDT >=90 and every neighboring residue pLDDT >=70 are required. A low target or neighbor excludes that option from this local gate; missing/ambiguous evidence holds it.
5. Export neighbor >=90 as the previously frozen sensitivity flag only. Zero other-residue neighbors yield a blank minimum (not zero) and a vacuously satisfied neighbor rule; target and all remaining context requirements still apply.
6. Parse the linked PAE matrix for model-relative dimensions and finite nonnegative values, accounting for integer-rounded matrices (e.g. 32 with reported maximum 31.75). Export the maximum of both target-to-neighbor and neighbor-to-target entries. PAE is descriptive: no newly selected PAE cutoff or score is introduced. Unresolved PAE is explicit and never imputed. Its interpretation belongs to later context review.

## Remaining context and meaning of pass

local_gate_status=pass means only that modeled local confidence and coherent target atom names pass these checks. structural_eligibility remains held for every mapped option. Full canonical versus fragment coverage is explicit; fragment boundaries are flagged. Neighbor atom completeness, native state, interdomain placement, partner/ligand/assembly context, experimental validation and chemical state are not approved here. No best model, biological assembly or feature-ready preparation is selected. No SASA, pKa or oxidation score is calculated. Missing/low-confidence structures do not demonstrate oxidative resistance.

Site aggregation is disjoint: a local-pass option is held pending context; unresolved local/experimental/mapping evidence is held; all mapped options failing the confidence gate with no remaining held evidence is excluded for this gate; no assessable mapping is unassessable. Experimental alternatives prevent declaring an entire site excluded just because its predicted option fails. Every original site/protein remains reportable.

## Outputs and errors

results/05b3_local_quality contains site_local_quality.csv (all 14,310 mapped options), predicted_model_audit.csv (assessed AlphaFold candidates), site_step_audit.csv (all 10,799 sites), protein_quality_coverage.csv (all 760 proteins), raw_file_checksums.csv, input_checksums.csv, summary.json, quality_report.md and SUCCESS.txt/FAILURE.txt. Local-pass option counts, unique site counts and protein counts are separate. Failure invalidates SUCCESS even if prior tables exist. SUCCESS means STEP05B3_GENERATED_REVIEW_PENDING; workstation outputs require review before the next binding.

## Validation

All 90 regression tests passed in the R 4.3.3 / DESeq2 1.42.0 regression runtime with Gemmi 0.7.5; the 16 local-quality tests were rerun after final validation refinements. Coverage includes exact confidence/distance boundaries, low target/neighbors, H/D exclusion, target-residue exclusion, coherent atom gaps, CIF/JSON conflicts, alternate/occupancy holds, fragment-relative numbering, PAE integer rounding, empty PAE without imputation, geometry/provenance mismatch rejection, frozen full input gates, experimental holds and full-universe audit retention.

A real accepted-hash AF-A0A087WRH0-F1 sample maps 35 sites: 21 pass the local gate and 14 fail target pLDDT. All neighborhood sets and gate decisions match an independent all-atom brute-force calculation. An end-to-end fixture processes that model's three raw files plus 11gl, retaining all 760 proteins / 10,799 sites; eight experimental mapped options remain held. Corrupting a raw byte correctly invalidates SUCCESS. These are representative validations, not the full 731-model workstation result.

## Primary documentation

- AlphaFold DB file formats, pLDDT/PAE distinction: https://github.com/google-deepmind/alphafold/blob/main/afdb/README.md
- AlphaFold confidence JSON serializer: https://github.com/google-deepmind/alphafold/blob/main/alphafold/common/confidence.py
- pLDDT and limits of interdomain confidence: https://www.ebi.ac.uk/training/online/courses/alphafold/inputs-and-outputs/evaluating-alphafolds-predicted-structures-using-confidence-scores/plddt-understanding-local-confidence/

## Revision 1.0.1 — preserve all upstream uncertainty in the site audit

Actual initial Ubuntu output at 7bb88bf0b17542d6218454c0c098f1b75fd476db passes option confidence/provenance checks: 4,580 predicted sites in 667 proteins pass the local gate; 5,036 options fail target confidence and 125 fail neighbor confidence. The initial site aggregator considered mapped options only and a prior reason code prioritizing mapped evidence. It could classify a site as excluded despite unresolved original experimental mapping or offered metadata.

The complete frozen mapped/unmapped candidate evidence now supplies unresolved_mapping_options and unresolved_metadata_review to site_step_audit.csv. A failed predicted option cannot exclude the whole site while such evidence is pending. On the uploaded data, 2,042 original excluded sites must become held; predicted local-gate decisions and the 4,580/667 counts do not change. Expected terminal counts after rerun: 6,964 held, 3,042 excluded, 793 unassessable. This is an audit-state correction, not a changed threshold or rescued/scored structure.

All 92 regression tests passed, including two specific regressions for failed predicted options with held experimental mappings or unresolved metadata. Initial actual output is preserved in docs/audits/ubuntu_local_quality_initial_2026-10-09 with correction-required review. Rerun the same Ubuntu command to generate consistent version 1.0.1 outputs before downstream binding.
