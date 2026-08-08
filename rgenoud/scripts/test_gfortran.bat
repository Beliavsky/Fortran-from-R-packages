@echo off
setlocal
cd /d "%~dp0\.."
set FFLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
if exist build_gfortran rmdir /s /q build_gfortran
mkdir build_gfortran

gfortran %FFLAGS% -J build_gfortran -I build_gfortran -c src\rgenoud_kinds.f90 -o build_gfortran\rgenoud_kinds.o || exit /b 1
gfortran %FFLAGS% -J build_gfortran -I build_gfortran -c src\rgenoud_random.f90 -o build_gfortran\rgenoud_random.o || exit /b 1
gfortran %FFLAGS% -J build_gfortran -I build_gfortran -c src\rgenoud_types.f90 -o build_gfortran\rgenoud_types.o || exit /b 1
gfortran %FFLAGS% -J build_gfortran -I build_gfortran -c src\rgenoud_derivatives.f90 -o build_gfortran\rgenoud_derivatives.o || exit /b 1
gfortran %FFLAGS% -J build_gfortran -I build_gfortran -c src\rgenoud_operators.f90 -o build_gfortran\rgenoud_operators.o || exit /b 1
gfortran %FFLAGS% -J build_gfortran -I build_gfortran -c src\rgenoud_core.f90 -o build_gfortran\rgenoud_core.o || exit /b 1
gfortran %FFLAGS% -J build_gfortran -I build_gfortran -c src\rgenoud_stats.f90 -o build_gfortran\rgenoud_stats.o || exit /b 1
gfortran %FFLAGS% -J build_gfortran -I build_gfortran -c src\rgenoud.f90 -o build_gfortran\rgenoud.o || exit /b 1

set OBJS=build_gfortran\rgenoud_kinds.o build_gfortran\rgenoud_random.o build_gfortran\rgenoud_types.o build_gfortran\rgenoud_derivatives.o build_gfortran\rgenoud_operators.o build_gfortran\rgenoud_core.o build_gfortran\rgenoud_stats.o build_gfortran\rgenoud.o
for %%F in (test\*.f90) do (
  gfortran %FFLAGS% -J build_gfortran -I build_gfortran "%%F" %OBJS% -o "build_gfortran\%%~nF.exe" || exit /b 1
  "build_gfortran\%%~nF.exe" || exit /b 1
)
echo All tests passed.
