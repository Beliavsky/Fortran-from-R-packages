@echo off
setlocal
cd /d %~dp0\..
if exist build-strict rmdir /s /q build-strict
mkdir build-strict
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
set MODS=-Jbuild-strict -Ibuild-strict

gfortran %FLAGS% %MODS% -c vendor\quadprog-fortran\src\quadprog_kinds.f90 -o build-strict\quadprog_kinds.o || exit /b 1
gfortran %FLAGS% %MODS% -c vendor\quadprog-fortran\src\quadprog_core.f90 -o build-strict\quadprog_core.o || exit /b 1
gfortran %FLAGS% %MODS% -c vendor\quadprog-fortran\src\quadprog.f90 -o build-strict\quadprog.o || exit /b 1
gfortran %FLAGS% %MODS% -c vendor\lpSolve-fortran\src\lpsolve_types.f90 -o build-strict\lpsolve_types.o || exit /b 1
gfortran %FLAGS% %MODS% -c vendor\lpSolve-fortran\src\lpsolve_simplex.f90 -o build-strict\lpsolve_simplex.o || exit /b 1
gfortran %FLAGS% %MODS% -c vendor\lpSolve-fortran\src\lpsolve_core.f90 -o build-strict\lpsolve_core.o || exit /b 1
gfortran %FLAGS% %MODS% -c vendor\lpSolve-fortran\src\lpsolve_special.f90 -o build-strict\lpsolve_special.o || exit /b 1
gfortran %FLAGS% %MODS% -c vendor\lpSolve-fortran\src\lpsolve.f90 -o build-strict\lpsolve.o || exit /b 1
gfortran %FLAGS% %MODS% -c src\limsolve_kinds.f90 -o build-strict\limsolve_kinds.o || exit /b 1
gfortran %FLAGS% %MODS% -c src\limsolve_types.f90 -o build-strict\limsolve_types.o || exit /b 1
gfortran %FLAGS% %MODS% -c src\limsolve_linalg.f90 -o build-strict\limsolve_linalg.o || exit /b 1
gfortran %FLAGS% %MODS% -c src\limsolve_inverse.f90 -o build-strict\limsolve_inverse.o || exit /b 1
gfortran %FLAGS% %MODS% -c src\limsolve_linear.f90 -o build-strict\limsolve_linear.o || exit /b 1
gfortran %FLAGS% %MODS% -c src\limsolve_ranges.f90 -o build-strict\limsolve_ranges.o || exit /b 1
gfortran %FLAGS% %MODS% -c src\limsolve_sampling.f90 -o build-strict\limsolve_sampling.o || exit /b 1
gfortran %FLAGS% %MODS% -c src\limsolve.f90 -o build-strict\limsolve.o || exit /b 1
set OBJS=build-strict\quadprog_kinds.o build-strict\quadprog_core.o build-strict\quadprog.o build-strict\lpsolve_types.o build-strict\lpsolve_simplex.o build-strict\lpsolve_core.o build-strict\lpsolve_special.o build-strict\lpsolve.o build-strict\limsolve_kinds.o build-strict\limsolve_types.o build-strict\limsolve_linalg.o build-strict\limsolve_inverse.o build-strict\limsolve_linear.o build-strict\limsolve_ranges.o build-strict\limsolve_sampling.o build-strict\limsolve.o
for %%F in (test\*.f90) do (
  gfortran %FLAGS% %MODS% %%F %OBJS% -o build-strict\%%~nF.exe || exit /b 1
  build-strict\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90) do (
  gfortran %FLAGS% %MODS% %%F %OBJS% -o build-strict\%%~nF.exe || exit /b 1
  build-strict\%%~nF.exe || exit /b 1
)
endlocal
