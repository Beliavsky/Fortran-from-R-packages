@echo off
rem SPDX-License-Identifier: GPL-3.0-only
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build-gfortran-windows
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\obj" "%BUILD%\mod" "%BUILD%\bin"
set FLAGS=-std=f2018 -Wall -Wextra -Werror -O2 -ffree-line-length-none -J"%BUILD%\mod" -I"%BUILD%\mod"
set OBJECTS=
for %%M in (bondvaluation_kinds bondvaluation_brazil_calendar bondvaluation_dates bondvaluation_daycount bondvaluation_schedule bondvaluation_pricing bondvaluation_compat bondvaluation) do (
  gfortran %FLAGS% -c "%ROOT%\src\%%M.f90" -o "%BUILD%\obj\%%M.o" || exit /b 1
  set OBJECTS=!OBJECTS! "%BUILD%\obj\%%M.o"
)
for %%F in ("%ROOT%\test\*.f90" "%ROOT%\app\*.f90" "%ROOT%\example\*.f90") do (
  gfortran %FLAGS% "%%F" !OBJECTS! -o "%BUILD%\bin\%%~nF.exe" || exit /b 1
)
for %%T in ("%BUILD%\bin\test_*.exe") do "%%T" || exit /b 1
"%BUILD%\bin\bondvaluation_demo.exe" >nul || exit /b 1
"%BUILD%\bin\day_count_comparison.exe" >nul || exit /b 1
"%BUILD%\bin\regular_bond.exe" >nul || exit /b 1
echo GNU Fortran Windows build: PASS
