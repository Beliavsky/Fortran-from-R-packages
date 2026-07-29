@echo off
setlocal EnableExtensions EnableDelayedExpansion
cd /d "%~dp0\.."

call :build build-validation-debug -O0
if errorlevel 1 exit /b 1
call :build build-validation-optimized -O2
if errorlevel 1 exit /b 1

echo validation: PASS
exit /b 0

:build
set "BDIR=%~1"
set "OPT=%~2"
if exist "%BDIR%" rmdir /s /q "%BDIR%"
mkdir "%BDIR%"

set "FLAGS=-std=f2018 %OPT% -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace -J %BDIR% -I %BDIR%"
set "MODULES=qrmtools_kinds qrmtools_types qrmtools_stats qrmtools_optimization qrmtools_distributions qrmtools_evt qrmtools_brownian qrmtools_black_scholes qrmtools_risk qrmtools_bounds qrmtools_returns qrmtools_hierarchy qrmtools_allocation qrmtools_tests qrmtools_garch qrmtools"
set "OBJECTS="

for %%M in (%MODULES%) do (
  gfortran %FLAGS% -c "src\%%M.f90" -o "%BDIR%\%%M.o"
  if errorlevel 1 exit /b 1
  set "OBJECTS=!OBJECTS! %BDIR%\%%M.o"
)

for %%F in (test\*.f90 app\*.f90 example\*.f90) do (
  gfortran %FLAGS% "%%F" !OBJECTS! -o "%BDIR%\%%~nF.exe"
  if errorlevel 1 exit /b 1
  "%BDIR%\%%~nF.exe"
  if errorlevel 1 exit /b 1
)
exit /b 0
