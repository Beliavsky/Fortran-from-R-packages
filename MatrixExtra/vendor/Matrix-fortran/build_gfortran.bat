@echo off
setlocal enabledelayedexpansion

set BUILD=build_gfortran
if exist %BUILD% rmdir /s /q %BUILD%
mkdir %BUILD%\mod
mkdir %BUILD%\bin
set FLAGS=-std=f2018 -O2 -Wall -Wextra -pedantic
set OBJECTS=

for %%F in (matrix_kinds matrix_status matrix_dense matrix_decompositions matrix_functions matrix_sparse matrix_sparse_solvers matrix_ordering matrix_io matrix_constructors matrix_sparse_stats matrix_advanced matrix) do (
  gfortran %FLAGS% -J %BUILD%\mod -I %BUILD%\mod -c src\%%F.f90 -o %BUILD%\%%F.o
  if errorlevel 1 exit /b 1
  set OBJECTS=!OBJECTS! %BUILD%\%%F.o
)

for %%F in (test\*.f90 app\*.f90 example\*.f90) do (
  gfortran %FLAGS% -I %BUILD%\mod %%F !OBJECTS! -o %BUILD%\bin\%%~nF.exe
  if errorlevel 1 exit /b 1
)

for %%F in (%BUILD%\bin\test_*.exe) do (
  %%F
  if errorlevel 1 exit /b 1
)

echo Build complete. Programs are in %BUILD%\bin.
