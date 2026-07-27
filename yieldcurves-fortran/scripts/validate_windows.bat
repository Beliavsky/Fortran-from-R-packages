@echo off
setlocal enabledelayedexpansion
cd /d %~dp0\..
if exist build\validation rmdir /s /q build\validation
mkdir build\validation\mod
mkdir build\validation\obj
mkdir build\validation\bin
set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Wconversion-extra -Werror -fcheck=all -fbacktrace -ffree-line-length-132
set OBJS=
for %%S in (yc_kinds yc_types yc_linalg yc_utils yc_splines yc_models yc_curve_ops yc_analysis yc_pca_mod yieldcurves) do (
  gfortran %FLAGS% -J build\validation\mod -I build\validation\mod -c src\%%S.f90 -o build\validation\obj\%%S.o
  if errorlevel 1 exit /b 1
  set OBJS=!OBJS! build\validation\obj\%%S.o
)
for %%T in (test\*.f90) do (
  gfortran %FLAGS% -J build\validation\mod -I build\validation\mod %%T !OBJS! -o build\validation\bin\%%~nT.exe
  if errorlevel 1 exit /b 1
  build\validation\bin\%%~nT.exe
  if errorlevel 1 exit /b 1
)
echo validation (Windows): PASS
