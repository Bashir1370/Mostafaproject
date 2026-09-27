param(
  [int]$Threads = 8
)

$ErrorActionPreference = "Stop"

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$rawDir = Join-Path $repoRoot "ferroptosis-only\data\raw\GSE247883"
$manifestFile = Join-Path $rawDir "ena_run_manifest.csv"
$refDir = Join-Path $repoRoot "ferroptosis-only\data\reference\ensembl101"
$indexDir = Join-Path $refDir "salmon_index"
$quantRoot = Join-Path $repoRoot "ferroptosis-only\data\processed\GSE247883\salmon"

New-Item -ItemType Directory -Force -Path $refDir | Out-Null
New-Item -ItemType Directory -Force -Path $quantRoot | Out-Null

if (-not (Test-Path $manifestFile)) {
  throw "Missing ENA manifest. Run 02a_GSE247883_prepare_raw.ps1 first."
}

if (-not (Get-Command salmon -ErrorAction SilentlyContinue)) {
  throw "salmon was not found on PATH. Install Salmon and reopen PowerShell."
}

if (-not (Get-Command curl.exe -ErrorAction SilentlyContinue)) {
  throw "curl.exe was not found."
}

$cdnaGz = Join-Path $refDir "Homo_sapiens.GRCh38.cdna.all.fa.gz"
$gtfGz = Join-Path $refDir "Homo_sapiens.GRCh38.101.gtf.gz"

$cdnaUrl = "https://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
$gtfUrl = "https://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.gtf.gz"

if (-not (Test-Path $cdnaGz)) {
  Write-Host "Downloading Ensembl release 101 cDNA reference..."
  & curl.exe -L --retry 5 -o $cdnaGz $cdnaUrl
  if ($LASTEXITCODE -ne 0) { throw "Failed to download cDNA reference." }
}

if (-not (Test-Path $gtfGz)) {
  Write-Host "Downloading Ensembl release 101 GTF..."
  & curl.exe -L --retry 5 -o $gtfGz $gtfUrl
  if ($LASTEXITCODE -ne 0) { throw "Failed to download GTF annotation." }
}

if (-not (Test-Path $indexDir)) {
  Write-Host "Building Salmon index..."
  & salmon index -t $cdnaGz -i $indexDir -p $Threads
  if ($LASTEXITCODE -ne 0) { throw "Salmon index failed." }
} else {
  Write-Host "Salmon index already exists; skipping."
}

$manifest = Import-Csv $manifestFile

foreach ($r in $manifest) {
  if (-not (Test-Path $r.fastq1_path) -or -not (Test-Path $r.fastq2_path)) {
    throw "Missing FASTQ files for $($r.gsm). Run 02a with -DownloadReads first."
  }

  $outDir = Join-Path $quantRoot $r.gsm
  $quantFile = Join-Path $outDir "quant.sf"

  if (Test-Path $quantFile) {
    Write-Host "$($r.gsm): quant.sf already exists; skipping."
    continue
  }

  Write-Host "Quantifying $($r.gsm) ..."
  & salmon quant -i $indexDir -l A -1 $r.fastq1_path -2 $r.fastq2_path --validateMappings --gcBias -p $Threads -o $outDir

  if ($LASTEXITCODE -ne 0) {
    throw "Salmon quantification failed for $($r.gsm)."
  }
}

Write-Host "All Salmon quantifications completed."
Write-Host "Next: run the R analysis script from repository root:"
Write-Host 'source("ferroptosis-only/scripts/02_GSE247883_rescue_replication.R")'
