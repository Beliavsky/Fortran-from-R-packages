@echo off
setlocal enabledelayedexpansion
if exist build\windows rmdir /s /q build\windows
mkdir build\windows\mod
mkdir build\windows\bin
set FLAGS=-std=f2018 -O0 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace
set SRC=src\nvmix_kinds.f90 src\nvmix_types.f90 src\nvmix_special.f90 src\nvmix_random.f90 src\nvmix_linalg.f90 src\nvmix_mixing.f90 src\nvmix_core.f90 src\nvmix_gamma_mix.f90 src\nvmix_distributions.f90 src\nvmix_risk_dependence.f90 src\nvmix_skewt.f90 src\nvmix_fitting.f90 src\nvmix_compat.f90 src\nvmix.f90
gfortran %FLAGS% -J build\windows\mod -I build\windows\mod -c %SRC% || exit /b 1
move *.o build\windows >nul
for %%F in (test\*.f90) do (
  gfortran %FLAGS% -J build\windows\mod -I build\windows\mod build\windows\*.o %%F -o build\windows\bin\%%~nF.exe || exit /b 1
  build\windows\bin\%%~nF.exe || exit /b 1
)
echo validation: PASS
