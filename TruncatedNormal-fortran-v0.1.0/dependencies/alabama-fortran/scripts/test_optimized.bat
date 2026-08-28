@echo off
setlocal
cd /d %~dp0\..
if exist build_optimized rmdir /s /q build_optimized
mkdir build_optimized\mod
set SRC=dependencies\numDeriv-fortran\src\numderiv_kinds.f90 dependencies\numDeriv-fortran\src\numderiv_types.f90 dependencies\numDeriv-fortran\src\numderiv_callbacks.f90 dependencies\numDeriv-fortran\src\numderiv_core.f90 dependencies\numDeriv-fortran\src\numderiv.f90 dependencies\roptim\src\roptim_lbfgsb_core.f90 dependencies\roptim\src\roptim.f90 src\alabama.f90
gfortran -std=f2018 -O3 -Wimplicit-interface -Werror=implicit-interface -J build_optimized\mod -I build_optimized\mod %SRC% test\test_alabama.f90 -o build_optimized\test_alabama.exe
if errorlevel 1 exit /b 1
build_optimized\test_alabama.exe
