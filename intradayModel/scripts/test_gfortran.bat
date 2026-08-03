@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build\gfortran-windows
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"

set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -fbacktrace -J. -I.
set SRC=%ROOT%\src
set TEST=%ROOT%\test

gfortran %FLAGS% -c ^
 "%SRC%\intraday_kinds.f90" ^
 "%SRC%\intraday_types.f90" ^
 "%SRC%\intraday_utils.f90" ^
 "%SRC%\intraday_kalman.f90" ^
 "%SRC%\intraday_fit.f90" ^
 "%SRC%\intraday_use.f90" ^
 "%SRC%\intraday_simulation.f90" ^
 "%SRC%\intraday_model.f90" ^
 "%TEST%\test_support.f90"
if errorlevel 1 exit /b 1

set OBJ=intraday_kinds.o intraday_types.o intraday_utils.o intraday_kalman.o intraday_fit.o intraday_use.o intraday_simulation.o intraday_model.o
for %%T in (test_kalman test_fit test_fixed_acceleration test_decompose) do (
  gfortran %FLAGS% %OBJ% test_support.o "%TEST%\%%T.f90" -o %%T.exe
  if errorlevel 1 exit /b 1
  %%T.exe
  if errorlevel 1 exit /b 1
)

echo All tests passed.
endlocal
