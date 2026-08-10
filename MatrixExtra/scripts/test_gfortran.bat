@echo off
setlocal enabledelayedexpansion
cd /d %~dp0\..
set FLAGS=-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
set BUILD=build_strict
if exist %BUILD% rmdir /s /q %BUILD%
mkdir %BUILD%\mod %BUILD%\obj %BUILD%\bin
set V=vendor\Matrix-fortran
for %%S in (matrix_kinds matrix_status matrix_dense matrix_decompositions matrix_functions matrix_sparse matrix_sparse_solvers matrix_ordering matrix_io matrix_constructors matrix_sparse_stats matrix_advanced matrix) do (
  gfortran %FLAGS% -J %BUILD%\mod -I %BUILD%\mod -c %V%\src\%%S.f90 -o %BUILD%\obj\%%S.o || exit /b 1
)
for %%S in (matrixextra_types matrixextra_conversions matrixextra_utils matrixextra_slice matrixextra_bind matrixextra_matmul matrixextra_ops matrixextra_linalg matrixextra_recycle matrixextra_pattern matrixextra) do (
  gfortran %FLAGS% -J %BUILD%\mod -I %BUILD%\mod -c src\%%S.f90 -o %BUILD%\obj\%%S.o || exit /b 1
)
set OBJS=
for %%O in (%BUILD%\obj\*.o) do set OBJS=!OBJS! %%O
for %%T in (test\*.f90) do (
  gfortran %FLAGS% -I %BUILD%\mod %%T !OBJS! -o %BUILD%\bin\%%~nT.exe || exit /b 1
  %BUILD%\bin\%%~nT.exe || exit /b 1
)
for %%T in (example\*.f90) do (
  gfortran %FLAGS% -I %BUILD%\mod %%T !OBJS! -o %BUILD%\bin\%%~nT.exe || exit /b 1
  %BUILD%\bin\%%~nT.exe || exit /b 1
)
endlocal
