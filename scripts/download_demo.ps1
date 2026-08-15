$ErrorActionPreference = 'Stop'
$destination = Join-Path $PSScriptRoot '..\data\raw\Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs.zip'
$destination = [System.IO.Path]::GetFullPath($destination)
$url = 'https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny/Xenium_Prime_MultiCellSeg_Mouse_Ileum_tiny_outs.zip'
New-Item -ItemType Directory -Force (Split-Path $destination) | Out-Null
Invoke-WebRequest -Uri $url -OutFile $destination
$hash = (Get-FileHash -Algorithm MD5 $destination).Hash.ToLower()
if ($hash -ne 'be9d917eaac2ade708c111132f0f379d') { throw "MD5 mismatch: $hash" }
Write-Host "Downloaded and verified: $destination"

