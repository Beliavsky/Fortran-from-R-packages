@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0\.."
if exist build_gfortran rmdir /s /q build_gfortran
mkdir build_gfortran
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -Jbuild_gfortran -Ibuild_gfortran
set SOURCES=src\fmbasics_kinds.f90 src\fmbasics_dates.f90 src\fmbasics_interpolation.f90 src\fmbasics_conventions.f90 src\fmbasics_rates.f90 src\fmbasics_curves.f90 src\fmbasics_credit.f90 src\fmbasics_volatility.f90 src\fmbasics_money.f90 src\fmbasics.f90
set OBJECTS=
for %%F in (%SOURCES%) do (
  gfortran %FLAGS% -c %%F -o build_gfortran\%%~nF.o || exit /b 1
  set OBJECTS=!OBJECTS! build_gfortran\%%~nF.o
)
for %%F in (test\*.f90) do (
  gfortran %FLAGS% %%F !OBJECTS! -o build_gfortran\%%~nF.exe || exit /b 1
  build_gfortran\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90 app\*.f90) do (
  gfortran %FLAGS% %%F !OBJECTS! -o build_gfortran\%%~nF.exe || exit /b 1
  build_gfortran\%%~nF.exe || exit /b 1
)
endlocal
