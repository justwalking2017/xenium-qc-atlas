$ErrorActionPreference = 'Stop'
$url = 'https://s3-us-west-2.amazonaws.com/10x.files/samples/xenium/3.0.0/Xenium_Prime_Breast_Cancer_FFPE/Xenium_Prime_Breast_Cancer_FFPE_outs.zip'
$destination = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..\data\raw\Xenium_Prime_Breast_Cancer_FFPE_outs.zip'))
$expectedBytes = 40984132655
New-Item -ItemType Directory -Force (Split-Path $destination) | Out-Null
curl.exe -L --fail --retry 8 --continue-at - --output $destination $url
$actualBytes = (Get-Item -LiteralPath $destination).Length
if ($actualBytes -ne $expectedBytes) { throw "Size mismatch: expected $expectedBytes, observed $actualBytes" }
Write-Host "Downloaded and size-verified: $destination"
