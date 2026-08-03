@echo off
setlocal EnableExtensions
set "ROOT=%~dp0.."
for %%I in ("%ROOT%") do set "ROOT=%%~fI"
set "BRIDGE=%ROOT%\rust_bridge"
set "BINDIR=%BRIDGE%\bin"
set "TARGETDIR=%BRIDGE%\target\release"

where cargo >nul 2>nul || (
  echo ERROR: Cargo/Rust is required to build the Clarabel backend.
  echo Install Rust from https://rustup.rs and rerun this script.
  exit /b 1
)
where tar >nul 2>nul || (
  echo ERROR: tar.exe is required to extract rust_bridge\vendor.tar.xz.
  echo Current Windows 10/11 installations normally include tar.exe.
  exit /b 1
)

if not exist "%BRIDGE%\vendor\clarabel" (
  echo Extracting vendored Rust dependencies...
  tar -xJf "%BRIDGE%\vendor.tar.xz" -C "%BRIDGE%" || exit /b 1
)
if not exist "%BRIDGE%\.cargo" mkdir "%BRIDGE%\.cargo"
(
  echo [source.crates-io]
  echo replace-with = "vendored-sources"
  echo [source.vendored-sources]
  echo directory = "vendor"
  echo [net]
  echo offline = true
) > "%BRIDGE%\.cargo\config.toml"

echo Building the Clarabel Rust bridge...
cargo build --manifest-path "%BRIDGE%\Cargo.toml" --release --offline --locked
if errorlevel 1 (
  echo The host Rust target failed. Trying x86_64-pc-windows-gnu...
  where rustup >nul 2>nul || (
    echo ERROR: rustup is needed to install the GNU Rust target.
    exit /b 1
  )
  rustup target add x86_64-pc-windows-gnu || exit /b 1
  cargo build --manifest-path "%BRIDGE%\Cargo.toml" --release --offline --locked --target x86_64-pc-windows-gnu || exit /b 1
  set "TARGETDIR=%BRIDGE%\target\x86_64-pc-windows-gnu\release"
)

if not exist "%TARGETDIR%\clarabel_fortran_bridge.dll" (
  echo ERROR: Cargo completed but clarabel_fortran_bridge.dll was not found.
  echo Expected: %TARGETDIR%\clarabel_fortran_bridge.dll
  exit /b 1
)
if not exist "%BINDIR%" mkdir "%BINDIR%"
copy /y "%TARGETDIR%\clarabel_fortran_bridge.dll" "%BINDIR%\clarabel_fortran_bridge.dll" >nul || exit /b 1

rem Rust's GNU Windows target can require MinGW runtime DLLs. Copy any that are
rem visible on PATH or bundled in the active Rust toolchain beside the bridge.
for %%D in (libgcc_s_seh-1.dll libwinpthread-1.dll libstdc++-6.dll libssp-0.dll) do (
  for /f "delims=" %%P in ('where %%D 2^>nul') do copy /y "%%P" "%BINDIR%\%%D" >nul
)
for /f "delims=" %%R in ('rustc --print sysroot 2^>nul') do set "RUST_SYSROOT=%%R"
if defined RUST_SYSROOT (
  for %%D in (libgcc_s_seh-1.dll libwinpthread-1.dll libstdc++-6.dll libssp-0.dll) do (
    for /r "%RUST_SYSROOT%" %%P in (%%D) do if not exist "%BINDIR%\%%D" copy /y "%%P" "%BINDIR%\%%D" >nul
  )
)

echo Built Clarabel backend:
echo   %BINDIR%\clarabel_fortran_bridge.dll

echo Verifying the backend through the Fortran interface...
where fpm >nul 2>nul
if errorlevel 1 (
  echo WARNING: fpm was not found, so the DLL could not be run-verified.
  echo Set CLARABEL_FORTRAN_BRIDGE to the DLL above when running your program.
) else (
  set "CLARABEL_FORTRAN_BRIDGE=%BINDIR%\clarabel_fortran_bridge.dll"
  cd /d "%ROOT%" || exit /b 1
  fpm run --target demo_clarabel
  if errorlevel 1 (
    echo ERROR: The bridge was built but Windows could not execute it.
    echo This usually means a dependent runtime DLL is missing.
    echo Inspect dependencies with: objdump -p "%BINDIR%\clarabel_fortran_bridge.dll" ^| findstr "DLL Name"
    exit /b 1
  )
)

echo Backend build and verification completed.
echo A normal fpm build no longer needs LIBRARY_PATH or -L flags.
endlocal
