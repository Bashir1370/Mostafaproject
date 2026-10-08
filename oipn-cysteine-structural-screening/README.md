# OIPN cysteine structural screening — protocol v0.1.0

Approved study design recorded on 2026-10-08. Target workstation: **Ubuntu Linux / Bash**.

> Among proteins corresponding to significantly differentially expressed protein-coding genes in mouse DRG after oxaliplatin, which contain cysteine sites with structural features consistent with susceptibility to thiol oxidation?

این پروژه پروتئین‌های متناظر با DEGs معنادار در DRG موش پس از اگزالی‌پلاتین را بر اساس جایگاه‌های سیستئینی اولویت‌بندی می‌کند. آستانهٔ اندازهٔ تغییر بیان نداریم. امتیاز ساختاری مستقل از log2FC و FDR است.

## Status

- Study scope, ten-step workflow and initial scoring rule: **recorded and frozen as v0.1.0**.
- Primary dataset: GSE286387 **candidate; current data/metadata audit pending**.
- Complementary dataset: GSE125002 **candidate; analyzed separately**.
- DEG analysis, structure retrieval and biological validation: **not run**.
- This commit contains a protocol, configuration, output contracts and Ubuntu instructions. It does not contain executable scientific analysis scripts or fabricated sample manifests/results.

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
