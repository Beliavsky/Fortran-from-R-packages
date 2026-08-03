@echo off
setlocal enabledelayedexpansion
cd /d %~dp0\..
if exist build-gfortran rmdir /s /q build-gfortran
mkdir build-gfortran\mod build-gfortran\obj build-gfortran\bin
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -O0 -g
for %%F in (nloptr_kinds nloptr_types nloptr_utils nloptr_derivatives nloptr_evaluation nloptr_solvers nloptr_api nloptr nloptr_example_functions) do (
  gfortran %FLAGS% -J build-gfortran\mod -I build-gfortran\mod -c src\%%F.f90 -o build-gfortran\obj\%%F.o || exit /b 1
)
for %%F in (test\*.f90 example\*.f90 app\*.f90) do (
  gfortran %FLAGS% -J build-gfortran\mod -I build-gfortran\mod %%F build-gfortran\obj\*.o -o build-gfortran\bin\%%~nF.exe || exit /b 1
  build-gfortran\bin\%%~nF.exe || exit /b 1
)
