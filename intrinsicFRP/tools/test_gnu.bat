@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0\.."
if exist build\windows rmdir /s /q build\windows
mkdir build\windows\mod
mkdir build\windows\obj
mkdir build\windows\bin
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -Wno-maybe-uninitialized -ffree-line-length-132 -O2
set SOURCES=intrinsicfrp_kinds intrinsicfrp_types intrinsicfrp_linalg intrinsicfrp_stats intrinsicfrp_hac intrinsicfrp_models intrinsicfrp_identification intrinsicfrp_oracle intrinsicfrp
for %%S in (%SOURCES%) do (
  gfortran %FLAGS% -Jbuild\windows\mod -Ibuild\windows\mod -c src\%%S.f90 -o build\windows\obj\%%S.o
  if errorlevel 1 exit /b 1
)
ar rcs build\windows\libintrinsicfrp.a build\windows\obj\*.o
if errorlevel 1 exit /b 1
for %%F in (test\*.f90 example\*.f90 app\*.f90) do (
  set NAME=%%~nF
  gfortran %FLAGS% -Ibuild\windows\mod %%F build\windows\libintrinsicfrp.a -o build\windows\bin\!NAME!.exe
  if errorlevel 1 exit /b 1
  build\windows\bin\!NAME!.exe
  if errorlevel 1 exit /b 1
)
echo GNU Fortran Windows validation: PASS
