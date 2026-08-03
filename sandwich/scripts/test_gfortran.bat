@echo off
setlocal enabledelayedexpansion
cd /d %~dp0\..
if exist build-gfortran rmdir /s /q build-gfortran
mkdir build-gfortran\mod build-gfortran\obj build-gfortran\bin
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -fbacktrace
set SOURCES=src\sandwich_kinds.f90 src\sandwich_status.f90 src\sandwich_utils.f90 src\sandwich_linalg.f90 src\sandwich_core.f90 src\sandwich_regression.f90 src\sandwich_kernels.f90 src\sandwich_auxiliary.f90 src\sandwich_hc.f90 src\sandwich_hac.f90 src\sandwich_cluster.f90 src\sandwich_panel.f90 src\sandwich_bootstrap.f90 src\sandwich.f90
for %%F in (%SOURCES%) do (
  gfortran %FLAGS% -Jbuild-gfortran\mod -Ibuild-gfortran\mod -c %%F -o build-gfortran\obj\%%~nF.o || exit /b 1
)
ar rcs build-gfortran\libsandwich.a build-gfortran\obj\*.o || exit /b 1
for %%F in (test\*.f90) do (
  gfortran %FLAGS% -Ibuild-gfortran\mod %%F build-gfortran\libsandwich.a -o build-gfortran\bin\%%~nF.exe || exit /b 1
  build-gfortran\bin\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90 app\*.f90) do (
  gfortran %FLAGS% -Ibuild-gfortran\mod %%F build-gfortran\libsandwich.a -o build-gfortran\bin\%%~nF.exe || exit /b 1
  build-gfortran\bin\%%~nF.exe >nul || exit /b 1
)
echo strict GNU Fortran validation: PASS
