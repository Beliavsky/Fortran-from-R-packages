@echo off
setlocal enabledelayedexpansion
cd /d %~dp0\..

set SOURCES=src\ghyp_kinds.f90 src\ghyp_special.f90 src\ghyp_linalg.f90 src\ghyp_rng.f90 src\ghyp_gig.f90 src\ghyp_model.f90 src\ghyp_distribution.f90 src\risksimul_types.f90 src\risksimul_math.f90 src\risksimul_portfolio.f90 src\risksimul_simulation.f90 src\risksimul.f90
set FLAGS=-std=f2018 -O0 -g -fcheck=all -fbacktrace -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror

if exist build_debug rmdir /s /q build_debug
mkdir build_debug
del /q *.o *.mod 2>nul

gfortran %FLAGS% -J build_debug -I build_debug -c %SOURCES% || exit /b 1
for %%F in (test\*.f90) do (
  gfortran %FLAGS% -J build_debug -I build_debug %%F *.o -o build_debug\%%~nF.exe || exit /b 1
  build_debug\%%~nF.exe || exit /b 1
)
for %%F in (app\*.f90 example\*.f90) do (
  gfortran %FLAGS% -J build_debug -I build_debug %%F *.o -o build_debug\%%~nF.exe || exit /b 1
  build_debug\%%~nF.exe >nul || exit /b 1
)

echo validation: PASS
