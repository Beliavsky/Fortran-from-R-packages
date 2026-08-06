@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0\.."
if exist build rmdir /s /q build
mkdir build\mod build\obj build\bin
set FLAGS=-std=f2018 -Wall -Wextra -Werror -O0 -g -fcheck=all -fbacktrace

gfortran %FLAGS% -Jbuild\mod -Ibuild\mod -c src\smith_wilson_kinds.f90 -o build\obj\smith_wilson_kinds.o || exit /b 1
gfortran %FLAGS% -Jbuild\mod -Ibuild\mod -c src\smith_wilson_linalg.f90 -o build\obj\smith_wilson_linalg.o || exit /b 1
gfortran %FLAGS% -Jbuild\mod -Ibuild\mod -c src\smith_wilson.f90 -o build\obj\smith_wilson.o || exit /b 1
gfortran %FLAGS% -Jbuild\mod -Ibuild\mod -c src\smith_wilson_yield_curve.f90 -o build\obj\smith_wilson_yield_curve.o || exit /b 1

for %%F in (test\test_*.f90) do (
  gfortran %FLAGS% -Ibuild\mod "%%F" build\obj\smith_wilson_kinds.o build\obj\smith_wilson_linalg.o build\obj\smith_wilson.o build\obj\smith_wilson_yield_curve.o -o "build\bin\%%~nF.exe" || exit /b 1
  "build\bin\%%~nF.exe" || exit /b 1
)
endlocal
