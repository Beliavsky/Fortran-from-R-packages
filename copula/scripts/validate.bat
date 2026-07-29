@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build-validation-windows
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"
set FLAGS=-std=f2018 -O0 -g -fcheck=all -fbacktrace -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror
set SOURCES=%ROOT%\src\copula_kinds.f90 %ROOT%\src\copula_types.f90 %ROOT%\src\copula_special.f90 %ROOT%\src\copula_random.f90 %ROOT%\src\copula_linalg.f90 %ROOT%\src\copula_families.f90 %ROOT%\src\copula_simulation.f90 %ROOT%\src\copula_dependence.f90 %ROOT%\src\copula_empirical.f90 %ROOT%\src\copula_fitting.f90 %ROOT%\src\copula_compositions.f90 %ROOT%\src\copula_special_discrete.f90 %ROOT%\src\copula.f90
gfortran %FLAGS% -J . -I . -c %SOURCES%
if errorlevel 1 exit /b 1
for %%F in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% -J . -I . "%%F" *.o -o "%%~nF.exe"
  if errorlevel 1 exit /b 1
  "%%~nF.exe"
  if errorlevel 1 exit /b 1
)
echo validation: PASS
