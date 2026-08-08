@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0\.."
if exist build-strict rmdir /s /q build-strict
mkdir build-strict\mod build-strict\obj build-strict\bin
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -Jbuild-strict\mod -Ibuild-strict\mod
for %%S in (gradient_kinds gradient_rng gradient_types gradient_stats gradient_sqgde gradient_benchmarks gradient) do (
  gfortran %FLAGS% -c src\%%S.f90 -o build-strict\obj\%%S.o || exit /b 1
)
for %%T in (test\*.f90) do (
  gfortran %FLAGS% %%T build-strict\obj\*.o -o build-strict\bin\%%~nT.exe || exit /b 1
  build-strict\bin\%%~nT.exe || exit /b 1
)
for %%E in (example\*.f90) do (
  gfortran %FLAGS% %%E build-strict\obj\*.o -o build-strict\bin\%%~nE.exe || exit /b 1
  build-strict\bin\%%~nE.exe || exit /b 1
)
endlocal
