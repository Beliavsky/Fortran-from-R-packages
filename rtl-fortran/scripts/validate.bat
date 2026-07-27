@echo off
rem SPDX-License-Identifier: MIT
setlocal enabledelayedexpansion
cd /d "%~dp0\.."

call :build_and_run build_validate_debug "-O0 -fcheck=all"
if errorlevel 1 exit /b 1
call :build_and_run build_validate_opt "-O2"
if errorlevel 1 exit /b 1

echo validation: PASS
exit /b 0

:build_and_run
set BUILD=%~1
set OPT=%~2
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod"
mkdir "%BUILD%\obj"
set FLAGS=-std=f2018 %OPT% -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -fbacktrace

for %%U in (rtl_kinds rtl_types rtl_stats rtl_calendar rtl_options rtl_processes rtl_fixed_income rtl_portfolio rtl_market rtl) do (
  gfortran %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" -c "src\%%U.f90" -o "%BUILD%\obj\%%U.o"
  if errorlevel 1 exit /b 1
)

set OBJECTS=
for %%U in (rtl_kinds rtl_types rtl_stats rtl_calendar rtl_options rtl_processes rtl_fixed_income rtl_portfolio rtl_market rtl) do (
  set OBJECTS=!OBJECTS! "%BUILD%\obj\%%U.o"
)

for %%T in (test\test_*.f90) do (
  gfortran %FLAGS% -I"%BUILD%\mod" "%%T" !OBJECTS! -o "%BUILD%\%%~nT.exe"
  if errorlevel 1 exit /b 1
  "%BUILD%\%%~nT.exe"
  if errorlevel 1 exit /b 1
)

for %%A in (app\*.f90 example\*.f90) do (
  gfortran %FLAGS% -I"%BUILD%\mod" "%%A" !OBJECTS! -o "%BUILD%\%%~nA.exe"
  if errorlevel 1 exit /b 1
  "%BUILD%\%%~nA.exe" >nul
  if errorlevel 1 exit /b 1
)
exit /b 0
