@echo off
setlocal EnableExtensions EnableDelayedExpansion

set MODE=%~1
if "%MODE%"=="" set MODE=checked
if "%FC%"=="" set FC=gfortran
set ROOT=%~dp0..
set BUILD=%ROOT%\build\%MODE%

if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\obj" || exit /b 1
mkdir "%BUILD%\bin" || exit /b 1

set COMMON=-std=f2018 -Wall -Wextra -Werror -ffree-line-length-none -J"%BUILD%\obj" -I"%BUILD%\obj"
if /i "%MODE%"=="checked" (
  set FLAGS=%COMMON% -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow
) else if /i "%MODE%"=="release" (
  set FLAGS=%COMMON% -O3
) else (
  echo Usage: %~nx0 [checked^|release]
  exit /b 2
)

set SOURCES=ptr_kinds ptr_types ptr_utils ptr_data ptr_indicators ptr_filters ptr_weighting ptr_performance ptr_backtest ptr_cross_sectional ptr_ml ptr_strategies ptr_optimization ptr_walk_forward portfolio_tester
set OBJECTS=
for %%S in (%SOURCES%) do (
  "%FC%" %FLAGS% -c "%ROOT%\src\%%S.f90" -o "%BUILD%\obj\%%S.o" || exit /b 1
  set OBJECTS=!OBJECTS! "%BUILD%\obj\%%S.o"
)

call :run_group test || exit /b 1
call :run_group example || exit /b 1
call :run_group app || exit /b 1
echo PortfolioTesteR %MODE% build: PASS
exit /b 0

:run_group
for %%F in ("%ROOT%\%~1\*.f90") do (
  "%FC%" %FLAGS% "%%~fF" %OBJECTS% -o "%BUILD%\bin\%%~nF.exe" || exit /b 1
  "%BUILD%\bin\%%~nF.exe" || exit /b 1
)
exit /b 0
