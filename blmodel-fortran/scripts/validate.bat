@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0\.."

if exist build-validation rmdir /s /q build-validation
mkdir build-validation\mod
mkdir build-validation\obj
mkdir build-validation\bin

set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace -O0 -g -Jbuild-validation\mod -Ibuild-validation\mod
set OBJECTS=

for %%F in (src\blmodel_kinds.f90 src\blmodel_types.f90 src\blmodel_linalg.f90 src\blmodel_utils.f90 src\blmodel_distributions.f90 src\blmodel_equilibrium.f90 src\blmodel_posterior.f90 src\blmodel.f90) do (
  gfortran %FLAGS% -c %%F -o build-validation\obj\%%~nF.o || exit /b 1
  set OBJECTS=!OBJECTS! build-validation\obj\%%~nF.o
)

gfortran %FLAGS% -c test\test_support.f90 -o build-validation\obj\test_support.o || exit /b 1

for %%F in (test_distributions test_equilibrium test_posterior test_blmodel) do (
  gfortran %FLAGS% test\%%F.f90 build-validation\obj\test_support.o !OBJECTS! -o build-validation\bin\%%F.exe || exit /b 1
  build-validation\bin\%%F.exe || exit /b 1
)

for %%F in (blmodel_demo) do (
  gfortran %FLAGS% app\%%F.f90 !OBJECTS! -o build-validation\bin\%%F.exe || exit /b 1
  build-validation\bin\%%F.exe > nul || exit /b 1
)

for %%F in (basic_black_litterman view_distributions) do (
  gfortran %FLAGS% example\%%F.f90 !OBJECTS! -o build-validation\bin\%%F.exe || exit /b 1
  build-validation\bin\%%F.exe > nul || exit /b 1
)

echo validation: PASS
endlocal
