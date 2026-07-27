@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build-validation-windows
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"

set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Werror -fcheck=all -fbacktrace
set SRC=%ROOT%\src\bidask_kinds.f90 %ROOT%\src\bidask_types.f90 %ROOT%\src\bidask_statistics.f90 %ROOT%\src\bidask_rng.f90 %ROOT%\src\bidask_estimators.f90 %ROOT%\src\bidask_windows.f90 %ROOT%\src\bidask_simulation.f90 %ROOT%\src\bidask.f90

gfortran %FLAGS% -J . -I . -c %SRC% || exit /b 1
set OBJ=bidask_kinds.o bidask_types.o bidask_statistics.o bidask_rng.o bidask_estimators.o bidask_windows.o bidask_simulation.o bidask.o
for %%F in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% -J . -I . "%%F" %OBJ% -o "%%~nF.exe" || exit /b 1
  "%%~nF.exe" || exit /b 1
)
for %%F in ("%ROOT%\app\*.f90" "%ROOT%\example\*.f90") do (
  gfortran %FLAGS% -J . -I . "%%F" %OBJ% -o "%%~nF.exe" || exit /b 1
  "%%~nF.exe" >nul || exit /b 1
)
echo validation: PASS
