$ErrorActionPreference = "Stop"
$Root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Pkg = Join-Path $Root "original\scip-master"
$ScipSrc = Join-Path $Pkg "inst\scip"
$SoplexSrc = Join-Path $Pkg "inst\soplex"
$Config = Join-Path $Pkg "inst\config"
$Build = Join-Path $Root "vendor\build"
$Lib = Join-Path $Root "vendor\lib"
$Jobs = if ($env:JOBS) { $env:JOBS } else { "2" }
New-Item -ItemType Directory -Force -Path $Build,$Lib | Out-Null

$Created = @()
function Add-Stub([string]$Dir) {
    if (-not (Test-Path $Dir)) {
        New-Item -ItemType Directory -Force -Path $Dir | Out-Null
        Set-Content -Path (Join-Path $Dir "CMakeLists.txt") -Value "# generated build stub"
        $script:Created += $Dir
    }
}
Add-Stub (Join-Path $SoplexSrc "check")
foreach ($d in @("check","tests","doc","examples","applications")) { Add-Stub (Join-Path $ScipSrc $d) }

try {
    $env:CC = if ($env:CC) { $env:CC } else { "gcc" }
    $env:CXX = if ($env:CXX) { $env:CXX } else { "g++" }
    cmake -G "MinGW Makefiles" -S $SoplexSrc -B (Join-Path $Build "soplex") `
        -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=ON `
        -DBUILD_SHARED_LIBS=OFF -DZLIB=OFF -DGMP=OFF -DMPFR=OFF `
        -DBOOST=OFF -DPAPILO=OFF -DQUADMATH=OFF "-DCMAKE_CXX_FLAGS=-I$Config"
    cmake --build (Join-Path $Build "soplex") --target libsoplex -j $Jobs

    cmake -G "MinGW Makefiles" -S $ScipSrc -B (Join-Path $Build "scip") `
        -DCMAKE_BUILD_TYPE=Release -DCMAKE_POSITION_INDEPENDENT_CODE=ON `
        -DBUILD_SHARED_LIBS=OFF "-DSOPLEX_DIR=$(Join-Path $Build 'soplex')" `
        -DZLIB=OFF -DGMP=OFF -DREADLINE=OFF -DPAPILO=OFF -DZIMPL=OFF `
        -DAMPL=OFF -DIPOPT=OFF -DWORHP=OFF -DCONOPT=OFF -DLAPACK=OFF `
        -DAUTOBUILD=OFF -DSHARED=OFF -DEXPRINT=none -DLPS=spx -DTPI=none `
        "-DCMAKE_C_FLAGS=-I$Config -DSCIP_LONGINT_FORMAT=\`"lld\`"" `
        "-DCMAKE_CXX_FLAGS=-I$Config -DSCIP_LONGINT_FORMAT=\`"lld\`""
    cmake --build (Join-Path $Build "scip") --target libscip -j $Jobs

    & $env:CC -c (Join-Path $Root "csrc\scip_fortran_shim.c") `
        -o (Join-Path $Build "scip_fortran_shim.o") -std=c11 -O2 `
        "-I$ScipSrc\src" "-I$Build\scip" "-I$Config"
    & $env:CXX -c (Join-Path $Root "csrc\standalone_streams.cpp") `
        -o (Join-Path $Build "standalone_streams.o") -std=c++17 -O2 "-I$Config"

    $Mri = Join-Path $Build "merge.mri"
    @"
create $Lib\libscipfortran_backend.a
addlib $Build\scip\lib\libscip.a
addlib $Build\soplex\lib\libsoplex.a
addmod $Build\scip_fortran_shim.o
addmod $Build\standalone_streams.o
save
end
"@ | Set-Content $Mri
    cmd /c "ar -M < `"$Mri`""
    Write-Host "Built $Lib\libscipfortran_backend.a"
}
finally {
    foreach ($d in $Created) { Remove-Item -Recurse -Force $d -ErrorAction SilentlyContinue }
}
