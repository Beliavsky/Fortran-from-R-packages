@echo off
setlocal
cd /d %~dp0\..
where gfortran >nul 2>nul || (echo gfortran was not found & exit /b 1)
if exist build rmdir /s /q build
mkdir build
set SRC=src\streg_kinds.f90 src\streg_linalg.f90 src\streg_probability.f90 src\streg_optimize.f90 src\streg_core.f90 src\streg.f90
for %%D in (block_top conditional static dynamic inference_errors) do (
  echo == test\%%D ==
  gfortran -std=f2018 -ffree-line-length-none -Wall -Wextra -Wpedantic -Werror -O0 -g -fcheck=all -fbacktrace -J build -I build %SRC% test\%%D\test_support.f90 test\%%D\main.f90 -o build\test_%%D.exe || exit /b 1
  build\test_%%D.exe || exit /b 1
)
echo All checked tests passed.
endlocal
