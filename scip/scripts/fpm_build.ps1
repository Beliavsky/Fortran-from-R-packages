$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Backend = Join-Path $Root "vendor\lib\libscipfortran_backend.a"
if (-not (Test-Path $Backend)) { & (Join-Path $PSScriptRoot "build_vendor.ps1") }
$old = $env:LIBRARY_PATH
$env:LIBRARY_PATH = (Join-Path $Root "vendor\lib") + $(if ($old) { ";$old" } else { "" })
Push-Location $Root
try { fpm build @args } finally { Pop-Location }
