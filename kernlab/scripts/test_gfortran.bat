@echo off
setlocal
cd /d %~dp0\..
if exist build\gfortran-debug rmdir /s /q build\gfortran-debug
mkdir build\gfortran-debug\mod
mkdir build\gfortran-debug\bin
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow
set SRC=src\kernlab_kinds.f90 src\kernlab_types.f90 src\kernlab_linalg.f90 src\kernlab_kernels.f90 src\kernlab_core.f90 src\kernlab_unsupervised.f90 src\kernlab_mmd.f90 src\kernlab_supervised.f90 src\kernlab.f90
gfortran %FLAGS% -J build\gfortran-debug\mod -I build\gfortran-debug\mod -c %SRC%
move /y *.o build\gfortran-debug\ >nul
for %%F in (test\*.f90 example\*.f90 app\*.f90) do (
  gfortran %FLAGS% -J build\gfortran-debug\mod -I build\gfortran-debug\mod %%F build\gfortran-debug\*.o -o build\gfortran-debug\bin\%%~nF.exe
  if errorlevel 1 exit /b 1
  build\gfortran-debug\bin\%%~nF.exe
  if errorlevel 1 exit /b 1
)
endlocal
