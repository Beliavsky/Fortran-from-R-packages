@echo off
setlocal enabledelayedexpansion
cd /d %~dp0\..
if exist build-gfortran rmdir /s /q build-gfortran
mkdir build-gfortran\mod build-gfortran\obj build-gfortran\bin
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -Wimplicit-interface -fcheck=all -fbacktrace -O0
set SOURCES=src\maxlik_kinds.f90 src\maxlik_status.f90 src\maxlik_types.f90 src\maxlik_linalg.f90 src\maxlik_random.f90 src\maxlik_evaluation.f90 src\maxlik_solvers.f90 src\maxlik_inference.f90 src\maxlik_utilities.f90 src\maxlik_api.f90 src\maxlik.f90
for %%F in (%SOURCES%) do (
  gfortran %FLAGS% -Jbuild-gfortran\mod -Ibuild-gfortran\mod -c %%F -o build-gfortran\obj\%%~nF.o || exit /b 1
)
ar rcs build-gfortran\libmaxlik.a build-gfortran\obj\*.o || exit /b 1
for %%F in (test\*.f90) do (
  gfortran %FLAGS% -Ibuild-gfortran\mod %%F build-gfortran\libmaxlik.a -o build-gfortran\bin\%%~nF.exe || exit /b 1
  build-gfortran\bin\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90 app\*.f90) do (
  gfortran %FLAGS% -Ibuild-gfortran\mod %%F build-gfortran\libmaxlik.a -o build-gfortran\bin\%%~nF.exe || exit /b 1
  build-gfortran\bin\%%~nF.exe >nul || exit /b 1
)
echo strict GNU Fortran validation: PASS
