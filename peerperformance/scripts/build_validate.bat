@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0\.."
if exist build-validation rmdir /s /q build-validation
mkdir build-validation\mod build-validation\obj build-validation\bin
set FLAGS=-O0 -g -std=f2018 -Wall -Wextra -Werror -Wconversion-extra -fcheck=all -fbacktrace
set MODS=peerperformance_kinds peerperformance_types peerperformance_math peerperformance_linalg peerperformance_stats peerperformance_pi peerperformance_bootstrap peerperformance_screening peerperformance
set OBJS=
for %%M in (%MODS%) do (
  gfortran %FLAGS% -Jbuild-validation\mod -Ibuild-validation\mod -c src\%%M.f90 -o build-validation\obj\%%M.o || exit /b 1
  set OBJS=!OBJS! build-validation\obj\%%M.o
)
for %%T in (test\*.f90) do (
  gfortran %FLAGS% -Ibuild-validation\mod %%T !OBJS! -o build-validation\bin\%%~nT.exe || exit /b 1
  build-validation\bin\%%~nT.exe || exit /b 1
)
for %%T in (app\*.f90 example\*.f90) do (
  gfortran %FLAGS% -Ibuild-validation\mod %%T !OBJS! -o build-validation\bin\%%~nT.exe || exit /b 1
  build-validation\bin\%%~nT.exe >nul || exit /b 1
)
echo validation: PASS
