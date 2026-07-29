@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build-validation
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"

set SRC=%ROOT%\src\vasicekfit_kinds.f90 %ROOT%\src\vasicekfit_normal.f90 %ROOT%\src\vasicekfit_linalg.f90 %ROOT%\src\vasicekfit_distribution.f90 %ROOT%\src\vasicekfit_model.f90 %ROOT%\src\vasicekfit_inference.f90 %ROOT%\src\vasicekfit.f90
set OBJ=vasicekfit_kinds.o vasicekfit_normal.o vasicekfit_linalg.o vasicekfit_distribution.o vasicekfit_model.o vasicekfit_inference.o vasicekfit.o
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wconversion-extra -Wimplicit-interface -fcheck=all -fbacktrace -O0 -J . -I .

gfortran %FLAGS% -c %SRC%
if errorlevel 1 exit /b 1

for %%F in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% "%%F" %OBJ% -o "%%~nF.exe"
  if errorlevel 1 exit /b 1
  "%%~nF.exe"
  if errorlevel 1 exit /b 1
)
for %%F in ("%ROOT%\app\*.f90" "%ROOT%\example\*.f90") do (
  gfortran %FLAGS% "%%F" %OBJ% -o "%%~nF.exe"
  if errorlevel 1 exit /b 1
  "%%~nF.exe"
  if errorlevel 1 exit /b 1
)

echo validation: PASS
