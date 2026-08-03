@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0\.."
if exist build-gfortran rmdir /s /q build-gfortran
mkdir build-gfortran\mod build-gfortran\obj build-gfortran\bin
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -fbacktrace
set SOURCES=src\invgamstochvol_kinds.f90 src\invgamstochvol_status.f90 src\invgamstochvol_rng.f90 src\invgamstochvol_special.f90 src\invgamstochvol_model.f90 src\invgamstochvol.f90
for %%F in (%SOURCES%) do (
  gfortran %FLAGS% -Jbuild-gfortran\mod -Ibuild-gfortran\mod -c %%F -o build-gfortran\obj\%%~nF.o || exit /b 1
)
ar rcs build-gfortran\libinvgamstochvol.a build-gfortran\obj\invgamstochvol_kinds.o build-gfortran\obj\invgamstochvol_status.o build-gfortran\obj\invgamstochvol_rng.o build-gfortran\obj\invgamstochvol_special.o build-gfortran\obj\invgamstochvol_model.o build-gfortran\obj\invgamstochvol.o || exit /b 1
for %%F in (test\*.f90) do (
  gfortran %FLAGS% -Ibuild-gfortran\mod %%F build-gfortran\libinvgamstochvol.a -o build-gfortran\bin\%%~nF.exe || exit /b 1
  build-gfortran\bin\%%~nF.exe || exit /b 1
)
echo strict GNU Fortran validation: PASS
endlocal
