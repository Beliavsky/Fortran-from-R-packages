@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build_strict
set FLAGS=-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\matrix\mod" "%BUILD%\piqp\mod" "%BUILD%\bin"
set MSRC=%ROOT%\vendor\Matrix-fortran\src
for %%S in (matrix_kinds matrix_status matrix_dense matrix_decompositions matrix_functions matrix_sparse matrix_sparse_solvers matrix_ordering matrix_io matrix_constructors matrix_sparse_stats matrix_advanced matrix) do (
  gfortran %FLAGS% -J"%BUILD%\matrix\mod" -I"%BUILD%\matrix\mod" -c "%MSRC%\%%S.f90" -o "%BUILD%\matrix\%%S.o" || exit /b 1
)
for %%S in (piqp_kinds piqp_types piqp_linalg piqp_solver) do (
  gfortran %FLAGS% -J"%BUILD%\piqp\mod" -I"%BUILD%\piqp\mod" -c "%ROOT%\src\%%S.f90" -o "%BUILD%\piqp\%%S.o" || exit /b 1
)
gfortran %FLAGS% -J"%BUILD%\piqp\mod" -I"%BUILD%\piqp\mod" -I"%BUILD%\matrix\mod" -c "%ROOT%\src\piqp_matrix_adapter.f90" -o "%BUILD%\piqp\piqp_matrix_adapter.o" || exit /b 1
gfortran %FLAGS% -J"%BUILD%\piqp\mod" -I"%BUILD%\piqp\mod" -I"%BUILD%\matrix\mod" -c "%ROOT%\src\piqp.f90" -o "%BUILD%\piqp\piqp.o" || exit /b 1
set OBJS=
for %%O in ("%BUILD%\piqp\*.o") do set OBJS=!OBJS! "%%O"
for %%O in ("%BUILD%\matrix\*.o") do set OBJS=!OBJS! "%%O"
for %%F in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% -I"%BUILD%\piqp\mod" -I"%BUILD%\matrix\mod" "%%F" !OBJS! -o "%BUILD%\bin\%%~nF.exe" || exit /b 1
  "%BUILD%\bin\%%~nF.exe" || exit /b 1
)
for %%F in ("%ROOT%\example\*.f90") do (
  gfortran %FLAGS% -I"%BUILD%\piqp\mod" -I"%BUILD%\matrix\mod" "%%F" !OBJS! -o "%BUILD%\bin\%%~nF.exe" || exit /b 1
  "%BUILD%\bin\%%~nF.exe" || exit /b 1
)
endlocal
