@echo off
setlocal enabledelayedexpansion
cd /d %~dp0\..
if exist build-gfortran rmdir /s /q build-gfortran
mkdir build-gfortran\mod build-gfortran\obj build-gfortran\bin
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -fbacktrace
set SOURCES=src\jumptest_kinds.f90 src\jumptest_status.f90 src\jumptest_rng.f90 src\jumptest_probability.f90 src\jumptest_statistics.f90 src\jumptest_simulation.f90 src\jumptest.f90
for %%S in (%SOURCES%) do (
  gfortran %FLAGS% -Jbuild-gfortran\mod -Ibuild-gfortran\mod -c %%S -o build-gfortran\obj\%%~nS.o || exit /b 1
)
ar rcs build-gfortran\libjumptest.a build-gfortran\obj\*.o || exit /b 1
for %%T in (test\*.f90) do (
  gfortran %FLAGS% -Ibuild-gfortran\mod %%T build-gfortran\libjumptest.a -o build-gfortran\bin\%%~nT.exe || exit /b 1
  build-gfortran\bin\%%~nT.exe || exit /b 1
)
for %%T in (example\*.f90 app\*.f90) do (
  gfortran %FLAGS% -Ibuild-gfortran\mod %%T build-gfortran\libjumptest.a -o build-gfortran\bin\%%~nT.exe || exit /b 1
  build-gfortran\bin\%%~nT.exe >nul || exit /b 1
)
echo strict GNU Fortran validation: PASS
