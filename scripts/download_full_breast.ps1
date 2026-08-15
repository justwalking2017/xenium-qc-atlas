$ErrorActionPreference = 'Stop'
$url = 'https://s3-us-west-2.amazonaws.com/10x.files/samples/xenium/1.0.2/Xenium_V1_FFPE_Human_Breast_IDC/Xenium_V1_FFPE_Human_Breast_IDC_outs.zip'
$destination = Join-Path $PSScriptRoot '..\data\raw\Xenium_V1_FFPE_Human_Breast_IDC_outs.zip'
$destination = [System.IO.Path]::GetFullPath($destination)
$expectedBytes = 24352808439
New-Item -ItemType Directory -Force (Split-Path $destination) | Out-Null
curl.exe -L --fail --retry 5 --continue-at - --output $destination $url
$actualBytes = (Get-Item -LiteralPath $destination).Length
if ($actualBytes -ne $expectedBytes) { throw "Size mismatch: expected $expectedBytes, observed $actualBytes" }
Write-Host "Downloaded and size-verified: $destination"
