@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build\gfortran-debug
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"
set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Werror -fcheck=all -fbacktrace
for %%S in (quadprog_kinds quadprog_core quadprog infoset_kinds infoset_status infoset_types infoset_stats infoset_mixture infoset_core infoset_portfolio infoset) do (
  gfortran %FLAGS% -J. -c "%ROOT%\src\%%S.f90" || exit /b 1
)
set OBJECTS=quadprog_kinds.o quadprog_core.o quadprog.o infoset_kinds.o infoset_status.o infoset_types.o infoset_stats.o infoset_mixture.o infoset_core.o infoset_portfolio.o infoset.o
for %%F in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% -I. "%%F" %OBJECTS% -o "%%~nF.exe" || exit /b 1
  "%%~nF.exe" || exit /b 1
)
for %%F in ("%ROOT%\example\*.f90" "%ROOT%\app\*.f90") do (
  gfortran %FLAGS% -I. "%%F" %OBJECTS% -o "%%~nF.exe" || exit /b 1
  "%%~nF.exe" || exit /b 1
)
endlocal
