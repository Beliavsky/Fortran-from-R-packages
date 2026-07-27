@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0\.."

if exist build-validation-windows rmdir /s /q build-validation-windows
mkdir build-validation-windows\mod
mkdir build-validation-windows\obj
mkdir build-validation-windows\bin

set FLAGS=-std=f2018 -O0 -g -fcheck=all -fbacktrace -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -Jbuild-validation-windows\mod -Ibuild-validation-windows\mod

for %%F in (epo_kinds epo_types epo_linalg epo_statistics epo_core epo) do (
  gfortran %FLAGS% -c src\%%F.f90 -o build-validation-windows\obj\%%F.o
  if errorlevel 1 exit /b 1
)

set OBJECTS=build-validation-windows\obj\epo_kinds.o build-validation-windows\obj\epo_types.o build-validation-windows\obj\epo_linalg.o build-validation-windows\obj\epo_statistics.o build-validation-windows\obj\epo_core.o build-validation-windows\obj\epo.o

for %%F in (test\*.f90 app\*.f90 example\*.f90) do (
  set NAME=%%~nF
  gfortran %FLAGS% "%%F" %OBJECTS% -o build-validation-windows\bin\!NAME!.exe
  if errorlevel 1 exit /b 1
  build-validation-windows\bin\!NAME!.exe
  if errorlevel 1 exit /b 1
)

endlocal
