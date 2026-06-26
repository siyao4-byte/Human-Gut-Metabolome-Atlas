$ErrorActionPreference = "Stop"

$root = (Resolve-Path ".").Path
$inputDir = Join-Path $root "outputs\figure4\library-inputs"
$libraryDir = Join-Path $root "outputs\figure4\libraries"
$reportDir = Join-Path $root "outputs\figure4\runtime-reports"
New-Item -ItemType Directory -Force -Path $libraryDir, $reportDir | Out-Null

$libraries = @(
    @{ platform = "PheHex"; mode = "pos"; ion = "P" },
    @{ platform = "PheHex"; mode = "neg"; ion = "N" },
    @{ platform = "HILIC"; mode = "pos"; ion = "P" },
    @{ platform = "HILIC"; mode = "neg"; ion = "N" },
    @{ platform = "SAX"; mode = "pos"; ion = "P" },
    @{ platform = "SAX"; mode = "neg"; ion = "N" }
)

foreach ($library in $libraries) {
    $stem = "$($library.platform)-$($library.mode)-paper2-rt075"
    python "scripts\07_figure4\02_build_msp_library.py" `
        --input (Join-Path $inputDir "$stem.csv") `
        --output (Join-Path $libraryDir "$stem.msp") `
        --ion-mode $library.ion `
        --report (Join-Path $reportDir "$stem.json")
    if ($LASTEXITCODE -ne 0) {
        throw "MSP generation failed for $stem"
    }
}

$reports = Get-ChildItem -LiteralPath $reportDir -Filter "*.json" |
    ForEach-Object { Get-Content -LiteralPath $_.FullName -Raw | ConvertFrom-Json }
$reports |
    Select-Object input, output, ion_mode, input_rows, unique_usable_usi,
        entries_written, failures, elapsed_seconds, elapsed_minutes |
    Export-Csv -NoTypeInformation -Path (Join-Path $root "outputs\tables\figure4-library-generation-runtime.csv")
$reports | Format-Table -AutoSize
