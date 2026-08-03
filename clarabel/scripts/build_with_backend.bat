@echo off
setlocal EnableExtensions
set "ROOT=%~dp0.."
for %%I in ("%ROOT%") do set "ROOT=%%~fI"
call "%ROOT%\scripts\build_backend.bat" || exit /b 1
set "CLARABEL_FORTRAN_BRIDGE=%ROOT%\rust_bridge\bin\clarabel_fortran_bridge.dll"
cd /d "%ROOT%" || exit /b 1
if "%~1"=="" (
  fpm build
) else (
  fpm %*
)
endlocal
