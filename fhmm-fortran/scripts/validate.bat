@echo off
setlocal enabledelayedexpansion
cd /d %~dp0\..
if exist build_validation rmdir /s /q build_validation
mkdir build_validation
set FLAGS=-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace
set SOURCES=src\fhmm_kinds.f90 src\fhmm_types.f90 src\fhmm_math.f90 src\fhmm_distributions.f90 src\fhmm_parameters.f90 src\fhmm_algorithms.f90 src\fhmm_hierarchical.f90 src\fhmm_optimize.f90 src\fhmm_estimation.f90 src\fhmm_diagnostics.f90 src\fhmm_calendar.f90 src\fhmm.f90
gfortran %FLAGS% -J build_validation -I build_validation -c %SOURCES% || exit /b 1
for %%F in (test\*.f90) do (
  set EXE=build_validation\%%~nF.exe
  gfortran %FLAGS% -J build_validation -I build_validation *.o %%F -o !EXE! || exit /b 1
  !EXE! || exit /b 1
)
for %%F in (app\*.f90 example\*.f90) do (
  set EXE=build_validation\%%~nF.exe
  gfortran %FLAGS% -J build_validation -I build_validation *.o %%F -o !EXE! || exit /b 1
  !EXE! > nul || exit /b 1
)
del /q *.o *.mod 2>nul
rmdir /s /q build_validation
echo validation: PASS
