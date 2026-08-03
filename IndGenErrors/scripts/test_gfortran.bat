@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0\.."
if exist build_gfortran rmdir /s /q build_gfortran
mkdir build_gfortran\mod build_gfortran\bin build_gfortran\obj
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow
for %%S in (indgen_kinds indgen_types indgen_special indgen_cvm_tables indgen_moebius indgen_core indgenerrors) do (
  gfortran %FLAGS% -J build_gfortran\mod -I build_gfortran\mod -c src\%%S.f90 -o build_gfortran\obj\%%S.o || exit /b 1
)
for %%T in (test\*.f90) do (
  gfortran %FLAGS% -I build_gfortran\mod build_gfortran\obj\*.o %%T -o build_gfortran\bin\%%~nT.exe || exit /b 1
  build_gfortran\bin\%%~nT.exe || exit /b 1
)
for %%T in (example\*.f90 app\*.f90) do (
  gfortran %FLAGS% -I build_gfortran\mod build_gfortran\obj\*.o %%T -o build_gfortran\bin\%%~nT.exe || exit /b 1
  build_gfortran\bin\%%~nT.exe || exit /b 1
)
endlocal
