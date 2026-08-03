@echo off
setlocal enabledelayedexpansion
cd /d %~dp0\..
if exist build-gfortran rmdir /s /q build-gfortran
mkdir build-gfortran\mod build-gfortran\obj build-gfortran\bin
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -fbacktrace
set SOURCES=src\moments_kinds.f90 src\moments_status.f90 src\moments_probability.f90 src\moments_statistics.f90 src\moments_transforms.f90 src\moments_tests.f90 src\moments.f90
for %%F in (%SOURCES%) do (
  gfortran %FLAGS% -Jbuild-gfortran\mod -Ibuild-gfortran\mod -c %%F -o build-gfortran\obj\%%~nF.o || exit /b 1
)
ar rcs build-gfortran\libmoments.a build-gfortran\obj\*.o || exit /b 1
for %%F in (test\*.f90) do (
  gfortran %FLAGS% -Ibuild-gfortran\mod %%F build-gfortran\libmoments.a -o build-gfortran\bin\%%~nF.exe || exit /b 1
  build-gfortran\bin\%%~nF.exe || exit /b 1
)
echo strict GNU Fortran validation: PASS
