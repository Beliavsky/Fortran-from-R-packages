@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build-gfortran
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -Wno-compare-reals -Wno-maybe-uninitialized -Wno-intrinsic-shadow -fcheck=all -ffpe-trap=invalid,zero,overflow -O0 -g
set MODS=pracma_kinds pracma_status pracma_callbacks pracma_types quadprog_kinds quadprog_core quadprog pracma_basic pracma_linalg pracma_polynomial pracma_special pracma_differentiation pracma_integration pracma_roots pracma_optimization pracma_interpolation pracma_ode pracma_signal_stats pracma_geometry pracma_combinatorics pracma_compat pracma
set OBJS=
for %%F in (%MODS%) do (
  gfortran %FLAGS% -c -J . -I . "%ROOT%\src\%%F.f90" || exit /b 1
  set OBJS=!OBJS! %%F.o
)
for %%T in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% -J . -I . "%%T" !OBJS! -o "%%~nT.exe" || exit /b 1
  "%%~nT.exe" || exit /b 1
)
endlocal
