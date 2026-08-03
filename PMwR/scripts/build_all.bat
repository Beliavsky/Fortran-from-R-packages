@echo off
setlocal enabledelayedexpansion
set MODE=%1
if "%MODE%"=="" set MODE=check
if "%FC%"=="" set FC=gfortran
set ROOT=%~dp0..
set BUILD=%ROOT%\build\%MODE%
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod"
mkdir "%BUILD%\bin"

set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic
if /I "%MODE%"=="check" set FLAGS=%FLAGS% -O0 -g -fcheck=all -fbacktrace
if /I "%MODE%"=="release" set FLAGS=%FLAGS% -O3
if /I not "%MODE%"=="check" if /I not "%MODE%"=="release" exit /b 2

set SOURCES=pmwr_kinds.f90 pmwr_types.f90 pmwr_utils.f90 pmwr_returns.f90 pmwr_portfolio.f90 pmwr_trades.f90 pmwr_analysis.f90 pmwr_backtest.f90 pmwr_identifiers.f90 pmwr.f90
set OBJECTS=
for %%S in (%SOURCES%) do (
  %FC% %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" -c "%ROOT%\src\%%S" -o "%BUILD%\%%~nS.o" || exit /b 1
  set OBJECTS=!OBJECTS! "%BUILD%\%%~nS.o"
)

for %%T in ("%ROOT%\test\*.f90") do (
  %FC% %FLAGS% -I"%BUILD%\mod" "%%T" !OBJECTS! -o "%BUILD%\bin\%%~nT.exe" || exit /b 1
  "%BUILD%\bin\%%~nT.exe" || exit /b 1
)
for %%E in ("%ROOT%\example\*.f90") do (
  %FC% %FLAGS% -I"%BUILD%\mod" "%%E" !OBJECTS! -o "%BUILD%\bin\%%~nE.exe" || exit /b 1
  "%BUILD%\bin\%%~nE.exe" >nul || exit /b 1
)
%FC% %FLAGS% -I"%BUILD%\mod" "%ROOT%\app\demo_pmwr.f90" !OBJECTS! -o "%BUILD%\bin\demo_pmwr.exe" || exit /b 1
"%BUILD%\bin\demo_pmwr.exe" >nul || exit /b 1

echo PMwR-fortran %MODE% build: PASS
endlocal
