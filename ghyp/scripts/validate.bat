@echo off
setlocal enabledelayedexpansion
cd /d %~dp0\..

set SRC=src\ghyp_kinds.f90 src\ghyp_special.f90 src\ghyp_rng.f90 src\ghyp_linalg.f90 src\ghyp_gig.f90 src\ghyp_model.f90 src\ghyp_distribution.f90 src\ghyp_risk.f90 src\ghyp_optimize.f90 src\ghyp_fitting.f90 src\ghyp_portfolio.f90 src\ghyp_utilities.f90 src\ghyp.f90

if exist build-validation-windows rmdir /s /q build-validation-windows
mkdir build-validation-windows
mkdir build-validation-windows\mod
mkdir build-validation-windows\obj
mkdir build-validation-windows\bin

gfortran -std=f2018 -O0 -Wall -Wextra -Wconversion-extra -Werror -fcheck=all -fbacktrace -J build-validation-windows\mod -I build-validation-windows\mod -c %SRC%
if errorlevel 1 exit /b 1
move /y *.o build-validation-windows\obj\ >nul

for %%F in (test\*.f90 app\*.f90 example\*.f90) do (
  set NAME=%%~nF
  gfortran -std=f2018 -O0 -Wall -Wextra -Wconversion-extra -Werror -fcheck=all -fbacktrace -I build-validation-windows\mod %%F build-validation-windows\obj\*.o -o build-validation-windows\bin\!NAME!.exe
  if errorlevel 1 exit /b 1
)

for %%F in (build-validation-windows\bin\test_*.exe) do (
  %%F
  if errorlevel 1 exit /b 1
)

build-validation-windows\bin\ghyp_demo.exe >nul
if errorlevel 1 exit /b 1
build-validation-windows\bin\distributions_and_risk.exe >nul
if errorlevel 1 exit /b 1
build-validation-windows\bin\fitting_and_portfolio.exe >nul
if errorlevel 1 exit /b 1

echo validation: PASS
endlocal
