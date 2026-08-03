@echo off
setlocal enabledelayedexpansion
cd /d %~dp0\..
if exist build-gfortran rmdir /s /q build-gfortran
mkdir build-gfortran\mod build-gfortran\obj build-gfortran\bin
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -fbacktrace
set SOURCES=src\gnorm_kinds.f90 src\gnorm_status.f90 src\gnorm_special.f90 src\gnorm_rng.f90 src\gnorm_distribution.f90 src\gnorm.f90
for %%S in (%SOURCES%) do (
  gfortran %FLAGS% -Jbuild-gfortran\mod -Ibuild-gfortran\mod -c %%S -o build-gfortran\obj\%%~nS.o || exit /b 1
)
ar rcs build-gfortran\libgnorm.a build-gfortran\obj\*.o || exit /b 1
for %%S in (test\*.f90) do (
  gfortran %FLAGS% -Ibuild-gfortran\mod %%S build-gfortran\libgnorm.a -o build-gfortran\bin\%%~nS.exe || exit /b 1
  build-gfortran\bin\%%~nS.exe || exit /b 1
)
for %%S in (example\*.f90 app\*.f90) do (
  gfortran %FLAGS% -Ibuild-gfortran\mod %%S build-gfortran\libgnorm.a -o build-gfortran\bin\%%~nS.exe || exit /b 1
  build-gfortran\bin\%%~nS.exe >nul || exit /b 1
)
echo strict GNU Fortran validation: PASS
