# Ubuntu Linux working guide

Target shell: Bash. Record the actual Ubuntu release, architecture, R, Python and package versions; no exact scientific dependency lock is claimed until the implementation is exercised. All new workflow commands will use Linux paths and Rscript/python3, not PowerShell.

## Retrieve repository

For a new checkout:

```bash
git clone https://github.com/Bashir1370/Mostafaproject.git
cd Mostafaproject
```

For an existing clean checkout, run `git pull --ff-only` from its root. Resolve any local changes before updating; do not discard them automatically.

## Basic system prerequisites

These commands are instructions for the user's Ubuntu machine; they were not executed there by this commit.

```bash
sudo apt update
sudo apt install git curl ca-certificates python3 python3-venv python3-pip r-base r-base-dev build-essential libcurl4-openssl-dev libssl-dev libxml2-dev
```

## Dedicated Python environment

From the repository root:

```bash
python3 -m venv oipn-cysteine-structural-screening/.venv
source oipn-cysteine-structural-screening/.venv/bin/activate
python -m pip install --upgrade pip
```

Scientific Python dependencies planned for later implementation: NumPy, pandas, Biopython, FreeSASA and PROPKA. Install/version-lock them when feature extraction is implemented and checked. A venv isolates Python dependencies; do not install into Ubuntu's system Python with sudo pip.

## Dedicated R library

```bash
mkdir -p oipn-cysteine-structural-screening/.r-library
export R_LIBS_USER="$PWD/oipn-cysteine-structural-screening/.r-library"
Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")'
Rscript -e 'BiocManager::install(c("DESeq2", "GEOquery", "AnnotationDbi", "org.Mm.eg.db"), ask=FALSE, update=FALSE)'
```

The final implementation will use one compatible R/Bioconductor release and record sessionInfo(). R packages needed only for later reporting/plots will be added when used, not treated as already tested dependencies. Annotation package symbols do not replace the release-specific gene/protein mapping audit.

## Check the recorded configuration

```bash
python3 -m json.tool oipn-cysteine-structural-screening/config/analysis_parameters.json
Rscript --version
python3 --version
lsb_release -ds
uname -m
```

The JSON command checks syntax only. The ten biological analysis scripts will be implemented stepwise and documented when ready; there is no currently runnable end-to-end analysis command.

## Directory and execution conventions

- Execute analysis from the repository root or use scripts resolving their own root; no dependence on a hard-coded personal home path.
- Config/docs are versioned. Runtime downloads, environments, caches and bulky results are excluded from Git within this workspace.
- Keep future data under data/raw, data/processed and data/structures; outputs under results/01 through results/10; logs under logs.
- Logs record input checksums, commands, software versions and protocol version. Downloads must be validated before use and resumed/retried without overwriting a valid file silently.
- No existing ferroptosis script is executed as part of this structural workflow. Its loaders are historical context only; this workflow needs its own audited mapping and outputs.

## Next executable task

Steps 01–02 are accepted and the user reports successful first-pass Step 03. Pull the latest commit and run the standard-library Python Step 03b resolver in the Ubuntu terminal; see [Step 03b](STEP_03B.md). It uses the completed Step 03 outputs/cache and HTTPS access to www.ebi.ac.uk. No extra Python/R package installation is needed.

## RStudio preparation verified by user (2026-10-08)

User workstation: R 4.3.3; Bioconductor 3.18; DESeq2 1.42.1. For Step 01b only DESeq2 and jsonlite are required. The QC script adds the dedicated project library to `.libPaths()` automatically, including when sourced from RStudio. It does not install or update packages during analysis.

To install the small JSON configuration reader in RStudio, if missing:

```r
project_lib <- "/home/bashir/Desktop/Mostafaproject/oipn-cysteine-structural-screening/.r-library"
dir.create(project_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(project_lib, .libPaths()))
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite", lib = project_lib, repos = "https://cloud.r-project.org")
}
source("/home/bashir/Desktop/Mostafaproject/oipn-cysteine-structural-screening/scripts/01b_qc_GSE286387.R")
```

Outputs are saved without changing the RStudio working directory. Provide `qc_report.md`, `sample_metrics.csv`, and `pca_top500.png` for review. Full diagnostics are collected in `QC_plots.pdf`. If FAILURE.txt exists or SUCCESS.txt is absent, the run is incomplete; rerun successfully before interpreting results. No DEG fitting occurs here.

## Step 04 offline inventory

After pulling the reviewed mapping snapshots and binding:

```bash
python3 oipn-cysteine-structural-screening/scripts/04_inventory_GSE286387.py
cat oipn-cysteine-structural-screening/results/04_cysteine_inventory/inventory_report.md
```

Python standard library suffices. No new R package, external download or structure tool is required. SUCCESS.txt must exist and FAILURE.txt must be absent. Send the report before proceeding to Step 05.
