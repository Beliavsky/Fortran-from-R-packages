@echo off
setlocal
cd /d "%~dp0\.."
if exist build-gfortran rmdir /s /q build-gfortran
mkdir build-gfortran\mod
mkdir build-gfortran\bin
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -fbacktrace
set SOURCES=src\fixedincome_kinds.f90 src\fixedincome_types.f90 src\fixedincome_terms.f90 src\fixedincome_compounding.f90 src\fixedincome_interpolation.f90 src\fixedincome_curves.f90 src\fixedincome.f90
for %%F in (test\*.f90) do (
  gfortran %FLAGS% -J build-gfortran\mod -I build-gfortran\mod %SOURCES% "%%F" -o "build-gfortran\bin\%%~nF.exe" || exit /b 1
  "build-gfortran\bin\%%~nF.exe" || exit /b 1
)
for %%F in (app\*.f90 example\*.f90) do (
  gfortran %FLAGS% -J build-gfortran\mod -I build-gfortran\mod %SOURCES% "%%F" -o "build-gfortran\bin\%%~nF.exe" || exit /b 1
  "build-gfortran\bin\%%~nF.exe" >nul || exit /b 1
)
echo strict GNU Fortran validation: PASS
