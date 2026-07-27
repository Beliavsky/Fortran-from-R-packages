@echo off
setlocal enabledelayedexpansion
if "%FC%"=="" set FC=gfortran
set FLAGS=-std=f2018 -Wall -Wextra -Wconversion-extra -Werror -Wno-compare-reals -O0 -fcheck=all -fbacktrace
if exist build-validation rmdir /s /q build-validation
mkdir build-validation
mkdir build-validation\mod
mkdir build-validation\obj
mkdir build-validation\bin

set SOURCES=src\highfrequency_kinds.f90 src\highfrequency_types.f90 src\highfrequency_stats.f90 src\highfrequency_linalg.f90 src\highfrequency_data.f90 src\highfrequency_cleaning.f90 src\highfrequency_realized.f90 src\highfrequency_optimize.f90 src\highfrequency_models.f90 src\highfrequency_jumps.f90 src\highfrequency_leadlag.f90 src\highfrequency_spot.f90 src\highfrequency_remedi.f90 src\highfrequency.f90
set OBJECTS=

for %%S in (%SOURCES%) do (
  set NAME=%%~nS
  %FC% %FLAGS% -J build-validation\mod -I build-validation\mod -c %%S -o build-validation\obj\!NAME!.o
  if errorlevel 1 exit /b 1
  set OBJECTS=!OBJECTS! build-validation\obj\!NAME!.o
)

for %%T in (test\*.f90) do (
  %FC% %FLAGS% -J build-validation\mod -I build-validation\mod %%T !OBJECTS! -o build-validation\bin\%%~nT.exe
  if errorlevel 1 exit /b 1
  build-validation\bin\%%~nT.exe
  if errorlevel 1 exit /b 1
)

for %%T in (app\*.f90 example\*.f90) do (
  %FC% %FLAGS% -J build-validation\mod -I build-validation\mod %%T !OBJECTS! -o build-validation\bin\%%~nT.exe
  if errorlevel 1 exit /b 1
  build-validation\bin\%%~nT.exe >nul
  if errorlevel 1 exit /b 1
)

echo validation: PASS
endlocal
