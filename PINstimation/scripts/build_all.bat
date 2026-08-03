@echo off
setlocal enabledelayedexpansion
set MODE=%1
if "%MODE%"=="" set MODE=check
set ROOT=%~dp0..
set BUILD=%ROOT%\build-%MODE%
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Wimplicit-interface -ffree-line-length-none
if /I "%MODE%"=="release" (
  set FLAGS=%FLAGS% -O3
) else (
  set FLAGS=%FLAGS% -O0 -g -fcheck=all -fbacktrace -finit-real=snan -finit-integer=-999999
)
for %%M in (pinstimation_kinds pinstimation_types pinstimation_math pinstimation_optimization pinstimation_pin pinstimation_mpin pinstimation_adjpin pinstimation_data pinstimation_vpin pinstimation) do (
  gfortran !FLAGS! -J"%BUILD%\mod" -I"%BUILD%\mod" -c "%ROOT%\src\%%M.f90" -o "%BUILD%\obj\%%M.o" || exit /b 1
)
set OBJECTS=
for %%O in ("%BUILD%\obj\*.o") do set OBJECTS=!OBJECTS! "%%O"
for %%S in ("%ROOT%\test\*.f90" "%ROOT%\example\*.f90" "%ROOT%\app\*.f90") do (
  gfortran !FLAGS! -I"%BUILD%\mod" "%%S" !OBJECTS! -o "%BUILD%\bin\%%~nS.exe" || exit /b 1
  "%BUILD%\bin\%%~nS.exe" || exit /b 1
)
echo All %MODE% builds and runs passed.
