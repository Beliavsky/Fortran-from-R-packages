@echo off
setlocal enabledelayedexpansion
cd /d %~dp0\..
if exist build_gfortran rmdir /s /q build_gfortran
mkdir build_gfortran\mod build_gfortran\bin build_gfortran\obj
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow
for %%F in (mixedind_kinds mixedind_types mixedind_rng mixedind_special mixedind_core mixedindtests) do (
  gfortran %FLAGS% -J build_gfortran\mod -I build_gfortran\mod -c src\%%F.f90 -o build_gfortran\obj\%%F.o || exit /b 1
)
for %%F in (test\*.f90) do (
  gfortran %FLAGS% -J build_gfortran\mod -I build_gfortran\mod build_gfortran\obj\*.o %%F -o build_gfortran\bin\%%~nF.exe || exit /b 1
  build_gfortran\bin\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90 app\*.f90) do (
  gfortran %FLAGS% -J build_gfortran\mod -I build_gfortran\mod build_gfortran\obj\*.o %%F -o build_gfortran\bin\%%~nF.exe || exit /b 1
  build_gfortran\bin\%%~nF.exe || exit /b 1
)
endlocal
