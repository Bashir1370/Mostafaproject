param(
  [switch]$DownloadReads
)

$ErrorActionPreference = "Stop"

$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$sampleFile = Join-Path $repoRoot "ferroptosis-only\config\GSE247883_samples.csv"
$outDir = Join-Path $repoRoot "ferroptosis-only\data\raw\GSE247883"
$fastqDir = Join-Path $outDir "fastq"

New-Item -ItemType Directory -Force -Path $outDir | Out-Null
New-Item -ItemType Directory -Force -Path $fastqDir | Out-Null

$samples = Import-Csv $sampleFile
$manifest = @()

Write-Host "Resolving ENA run/FASTQ metadata for GSE247883..."

foreach ($s in $samples) {
  $query = [uri]::EscapeDataString('experiment_accession="' + $s.srx + '"')
  $url = "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=$query&fields=run_accession,experiment_accession,fastq_ftp,fastq_md5,fastq_bytes,library_layout&format=tsv&limit=0"

  $txt = (Invoke-WebRequest -Uri $url -UseBasicParsing).Content
  $rows = $txt | ConvertFrom-Csv -Delimiter "`t"

  if ($rows.Count -ne 1) {
    throw "Expected exactly one ENA run for $($s.srx); found $($rows.Count)."
  }

  $r = $rows[0]

  if ($r.library_layout -ne "PAIRED") {
    throw "Expected PAIRED layout for $($s.srx), found '$($r.library_layout)'."
  }

  $ftp = $r.fastq_ftp -split ";"
  $md5 = $r.fastq_md5 -split ";"
  $bytes = $r.fastq_bytes -split ";"

  if ($ftp.Count -ne 2) {
    throw "Expected 2 FASTQ files for $($s.srx); found $($ftp.Count)."
  }

  $fq1Url = if ($ftp[0] -match '^https?://') { $ftp[0] } else { "https://$($ftp[0])" }
  $fq2Url = if ($ftp[1] -match '^https?://') { $ftp[1] } else { "https://$($ftp[1])" }

  $fq1 = Join-Path $fastqDir ([IO.Path]::GetFileName($ftp[0]))
  $fq2 = Join-Path $fastqDir ([IO.Path]::GetFileName($ftp[1]))

  $manifest += [pscustomobject]@{
    gsm = $s.gsm
    srx = $s.srx
    run = $r.run_accession
    condition = $s.condition
    replicate = $s.replicate
    fastq1_url = $fq1Url
    fastq2_url = $fq2Url
    fastq1_path = $fq1
    fastq2_path = $fq2
    fastq1_md5 = $md5[0]
    fastq2_md5 = $md5[1]
    fastq1_bytes = [double]$bytes[0]
    fastq2_bytes = [double]$bytes[1]
  }
}

$manifestFile = Join-Path $outDir "ena_run_manifest.csv"
$manifest | Export-Csv -NoTypeInformation -Path $manifestFile

$totalBytes = ($manifest | Measure-Object -Property fastq1_bytes -Sum).Sum +
              ($manifest | Measure-Object -Property fastq2_bytes -Sum).Sum
$totalGB = [math]::Round($totalBytes / 1GB, 2)

Write-Host ""
Write-Host "Manifest written to:"
Write-Host "  $manifestFile"
Write-Host "Total compressed FASTQ download size: $totalGB GB"
Write-Host ""

if (-not $DownloadReads) {
  Write-Host "Preflight only. No reads were downloaded."
  Write-Host "To download after reviewing the size:"
  Write-Host '  powershell -ExecutionPolicy Bypass -File "ferroptosis-only\scripts\02a_GSE247883_prepare_raw.ps1" -DownloadReads'
  exit 0
}

if (-not (Get-Command curl.exe -ErrorAction SilentlyContinue)) {
  throw "curl.exe was not found. It is normally included with modern Windows."
}

foreach ($r in $manifest) {
  foreach ($mate in 1,2) {
    $url = if ($mate -eq 1) { $r.fastq1_url } else { $r.fastq2_url }
    $path = if ($mate -eq 1) { $r.fastq1_path } else { $r.fastq2_path }
    $expectedMd5 = if ($mate -eq 1) { $r.fastq1_md5 } else { $r.fastq2_md5 }

    Write-Host "Downloading $($r.gsm) mate $mate ..."
    & curl.exe -L --retry 5 --retry-delay 5 -C - -o $path $url
    if ($LASTEXITCODE -ne 0) {
      throw "curl failed for $url"
    }

    $observedMd5 = (Get-FileHash -Algorithm MD5 -Path $path).Hash.ToLower()
    if ($observedMd5 -ne $expectedMd5.ToLower()) {
      throw "MD5 mismatch for $path"
    }
  }
}

Write-Host "All FASTQ files downloaded and MD5-validated."
