@echo off
rem SPDX-License-Identifier: Artistic-2.0
setlocal enabledelayedexpansion

set MODE=%1
if "%MODE%"=="" set MODE=strict
set ROOT=%~dp0..
set BUILD=%ROOT%\build\gfortran-%MODE%
set MODS=%BUILD%\mod
set OBJS=%BUILD%\obj
set BIN=%BUILD%\bin

if "%MODE%"=="strict" (
  set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Werror -Wno-compare-reals -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
) else if "%MODE%"=="release" (
  set FLAGS=-std=f2018 -O3 -Wall -Wextra -Werror -Wno-compare-reals
) else (
  echo Usage: %0 [strict^|release]
  exit /b 2
)

if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%MODS%" "%OBJS%" "%BIN%"

set SOURCES=ecd_kinds ecd_rng ecd_math ecd_core ecld_models ecd_processes ecd_options ecd_timeseries lamp_process ecd_fitting ecd_compat ecd_api
set OBJECTS=
for %%S in (%SOURCES%) do (
  gfortran %FLAGS% -J"%MODS%" -I"%MODS%" -c "%ROOT%\src\%%S.f90" -o "%OBJS%\%%S.o"
  if errorlevel 1 exit /b 1
  set OBJECTS=!OBJECTS! "%OBJS%\%%S.o"
)

for %%G in (test app example) do (
  for %%F in ("%ROOT%\%%G\*.f90") do (
    echo [%MODE%] %%G/%%~nF
    gfortran %FLAGS% -I"%MODS%" "%%F" !OBJECTS! -o "%BIN%\%%~nF.exe"
    if errorlevel 1 exit /b 1
    "%BIN%\%%~nF.exe"
    if errorlevel 1 exit /b 1
  )
)

echo All %MODE% builds and runs passed.
endlocal
