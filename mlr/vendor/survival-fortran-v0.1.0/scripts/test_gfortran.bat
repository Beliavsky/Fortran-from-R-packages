@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build-gfortran
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
for %%F in (splines_kinds splines_linalg splines_core splines_basis splines) do (
  gfortran %FLAGS% -c "%ROOT%\vendor\splines-fortran-v0.1.0\src\%%F.f90" || exit /b 1
)
for %%F in (survival_kinds survival_types survival_linalg survival_nonparametric survival_cox survival_aft survival_stats survival_utils survival_pspline survival) do (
  gfortran %FLAGS% -c "%ROOT%\src\%%F.f90" || exit /b 1
)
for %%T in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% "%%T" *.o -o test.exe || exit /b 1
  test.exe || exit /b 1
)
for %%E in ("%ROOT%\example\*.f90") do (
  gfortran %FLAGS% "%%E" *.o -o example.exe || exit /b 1
  example.exe || exit /b 1
)
endlocal
