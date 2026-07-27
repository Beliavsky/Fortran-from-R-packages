@echo off
setlocal enabledelayedexpansion
cd /d %~dp0\..
if exist build rmdir /s /q build
mkdir build\mod build\obj build\bin
set FLAGS=-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -O0 -fcheck=all -fbacktrace -Jbuild\mod -Ibuild\mod
set SOURCES=src\hdshop_kinds.f90 src\hdshop_linalg.f90 src\hdshop_stats.f90 src\hdshop_shrinkage.f90 src\hdshop_portfolio.f90 src\hdshop_inference.f90 src\hdshop_random.f90 src\hdshop_formulas.f90 src\hdshop_compat.f90 src\hdshop.f90
set OBJECTS=
for %%F in (%SOURCES%) do (
  set OBJ=build\obj\%%~nF.o
  gfortran %FLAGS% -c %%F -o !OBJ! || exit /b 1
  set OBJECTS=!OBJECTS! !OBJ!
)
for %%F in (test\*.f90) do (
  gfortran %FLAGS% %%F !OBJECTS! -o build\bin\%%~nF.exe || exit /b 1
  build\bin\%%~nF.exe || exit /b 1
)
for %%F in (app\*.f90 example\*.f90) do (
  gfortran %FLAGS% %%F !OBJECTS! -o build\bin\%%~nF.exe || exit /b 1
  build\bin\%%~nF.exe >nul || exit /b 1
)
echo validation: PASS
