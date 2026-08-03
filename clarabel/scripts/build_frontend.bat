@echo off
setlocal EnableExtensions EnableDelayedExpansion
set "ROOT=%~dp0.."
for %%I in ("%ROOT%") do set "ROOT=%%~fI"
set "BUILD=%ROOT%\build\frontend-windows"
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mock" || exit /b 1
cd /d "%BUILD%" || exit /b 1
set "FFLAGS=-std=f2018 -Wall -Wextra -Wpedantic -fcheck=all -fbacktrace -O0 -g"
set "CFLAGS=-std=c11 -Wall -Wextra -Wpedantic -O0 -g"

gcc !CFLAGS! -shared -I"%ROOT%\include" "%ROOT%\test\mock_backend.c" -o "%BUILD%\mock\clarabel_fortran_bridge.dll" || exit /b 1
gcc !CFLAGS! -c "%ROOT%\src\clarabel_dynamic_loader.c" || exit /b 1
for %%F in (clarabel_kinds.f90 clarabel_sparse.f90 clarabel_types.f90 clarabel_psd.f90 clarabel_c_api.f90 clarabel_solver.f90 clarabel.f90) do (
  gfortran !FFLAGS! -J. -I. -c "%ROOT%\src\%%F" || exit /b 1
)
set "CLARABEL_FORTRAN_BRIDGE=%BUILD%\mock\clarabel_fortran_bridge.dll"
for %%T in ("%ROOT%\test\*.f90") do (
  gfortran !FFLAGS! -J. -I. "%%T" *.o -o "%%~nT.exe" || exit /b 1
  "%%~nT.exe" || exit /b 1
)
for %%T in ("%ROOT%\example\*.f90" "%ROOT%\app\*.f90") do (
  gfortran !FFLAGS! -J. -I. "%%T" *.o -o "%%~nT.exe" || exit /b 1
)
unconstrained_qp.exe || exit /b 1
equality_qp.exe || exit /b 1
persistent_update.exe || exit /b 1
demo_clarabel.exe || exit /b 1
echo Frontend/runtime-loader checked validation passed.
endlocal
