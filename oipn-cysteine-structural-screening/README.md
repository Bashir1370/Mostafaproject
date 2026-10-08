# OIPN cysteine structural screening — protocol v0.1.0

Approved study design recorded on 2026-10-08. Target workstation: **Ubuntu Linux / Bash**.

> Among proteins corresponding to significantly differentially expressed protein-coding genes in mouse DRG after oxaliplatin, which contain cysteine sites with structural features consistent with susceptibility to thiol oxidation?

این پروژه پروتئین‌های متناظر با DEGs معنادار در DRG موش پس از اگزالی‌پلاتین را بر اساس جایگاه‌های سیستئینی اولویت‌بندی می‌کند. آستانهٔ اندازهٔ تغییر بیان نداریم. امتیاز ساختاری مستقل از log2FC و FDR است.

## Status

- Study scope, ten-step workflow and initial scoring rule: **recorded and frozen as v0.1.0**.
- Primary dataset: GSE286387 **Step 01 accepted after user QC and study-method review**; all ten samples retained with documented metadata limitations.
- Complementary dataset: GSE125002 **candidate; analyzed separately**.
- Step 02: **Ubuntu reproduction and uploaded snapshot reviewed**; 854 coding DEGs (477 up / 377 down), exact selected IDs match validation. [Review and limits](docs/audits/ubuntu_DE_2026-10-08/REVIEW.md); Step 03 input is checksum-bound in config/GSE286387_step03_input.json.
- Step 03: exact mouse gene mapping and canonical sequence retrieval implemented and tested; **workstation reproduction and held-gene review pending**.
- Structure retrieval and biological validation: **not run**.

## Ten steps

| Part | Steps | Purpose |
|---|---|---|
| Input | 01–02 | Dataset audit/QC and differential expression |
| Structural analysis | 03–09 | Mapping, cysteines, structures, chemical state, features, scores, validation |
| Interpretation | 10 | Biological interpretation and final candidate report |

## Documents

- [Detailed protocol](docs/PROTOCOL.md): every step, entry/exit criteria and outputs.
- [Scoring specification](docs/SCORING_SPEC.md): exact formula, ranks, ties and missing values.
- [Validation plan](docs/VALIDATION_PLAN.md): independent evidence and sensitivity analyses.
- [Ubuntu guide](docs/UBUNTU.md): local environment and execution conventions.
- [Output contracts](docs/OUTPUT_CONTRACTS.md): required tables and provenance.
- [Decision log](docs/DECISION_LOG.md): agreed choices and pending implementation decisions.
- [References](docs/REFERENCES.md): methodological and dataset sources.
- [Analysis parameters](config/analysis_parameters.json): machine-readable protocol settings.
- [Step registry](config/steps.csv): all ten steps.
- [Dataset candidates](config/dataset_candidates.csv): candidates, not an audited sample list.
- [Progress](STATUS.md): completed and next work.

## Ubuntu: retrieve this recorded design

From an existing checkout, run `git pull --ff-only`. For a new checkout:

```bash
git clone https://github.com/Bashir1370/Mostafaproject.git
cd Mostafaproject
python3 -m json.tool oipn-cysteine-structural-screening/config/analysis_parameters.json
```

The final command checks JSON syntax only; it does not run the study. See the Ubuntu guide before preparing R/Python dependencies.

## Interpretation boundaries

The score is a relative structural prioritization index, not a calibrated probability, oxidation rate, ROS concentration or proof of oxidation in OIPN. DEG restriction does not enumerate all oxidation-sensitive DRG proteins. RNA does not prove protein abundance or expressed isoform. High rank does not establish harmful consequences or causality in pain. This branch models free cysteine thiols; it does not cover every ROS target or chemical modification.

The workspace has its own configuration and outputs. It does not import ferroptosis signatures or existing classifier weights.

## Step 01: run on Ubuntu

```bash
cd /home/bashir/Desktop/Mostafaproject
git pull --ff-only
python3 oipn-cysteine-structural-screening/scripts/01_audit_GSE286387.py
cat oipn-cysteine-structural-screening/results/01_dataset_audit/qc_report.md
```

Python 3.8+ standard library only; no pip packages or R packages required for this download/integrity step. R sample QC is now implemented; DE testing follows only after QC/design review. See [Step 01 execution and findings](docs/STEP_01.md). Generated data/results are ignored by git. The committed audit snapshot is under `docs/audits/` and records the actual source checksums.

## Step 01b: R sample QC

After the audit completes, execute from the repository root:

```bash
Rscript oipn-cysteine-structural-screening/scripts/01b_qc_GSE286387.R
```

Or in the RStudio Console:

```r
source("/home/bashir/Desktop/Mostafaproject/oipn-cysteine-structural-screening/scripts/01b_qc_GSE286387.R")
```

The script loads the project's `.r-library` automatically. It requires DESeq2 and jsonlite (installation in [Ubuntu guide](docs/UBUNTU.md)). Outputs, including four PNGs and a combined PDF, are in `results/01b_sample_qc`. The QC script itself reports **QC_GENERATED_REVIEW_PENDING**; the subsequent review and accepted design are recorded separately in [design audit](docs/DATASET_DESIGN_AUDIT.md). See [QC specification and validation](docs/STEP_01B_QC.md).

## Step 02: differential expression

The accepted design is `~ condition`, oxaliplatin versus control, five libraries per group. BH uses all valid tested genes; significant coding genes proceed without a fold-change cutoff.

```bash
Rscript oipn-cysteine-structural-screening/scripts/02_DE_GSE286387.R
cat oipn-cysteine-structural-screening/results/02_differential_expression/de_report.md
```

Or in RStudio:

```r
source("/home/bashir/Desktop/Mostafaproject/oipn-cysteine-structural-screening/scripts/02_DE_GSE286387.R")
```

See [Step 02 inputs, outputs and validation](docs/STEP_02.md). A complete run requires SUCCESS.txt and absence of FAILURE.txt. Generated results remain ignored by Git; review and commit a deliberately selected snapshot after local reproduction.

## Step 03: canonical mouse protein mapping

```bash
python3 oipn-cysteine-structural-screening/scripts/03_map_GSE286387.py
cat oipn-cysteine-structural-screening/results/03_protein_mapping/mapping_report.md
```

Python 3.8+ standard library only. The checksum-bound Step 02 discovery snapshot is the input. Current exact Ensembl GeneId cross-references identify mouse UniProtKB entries; multiple accessions remain held instead of selecting an arbitrary reviewed/longest entry. Results include complete gene/candidate audits, unique-protein canonical FASTA and reference provenance. See [Step 03 criteria, cache and validation](docs/STEP_03.md). Large held-cohort coverage limitations must be reviewed before downstream inclusion is finalized.
