@echo off
setlocal
cd /d %~dp0\..
if exist build_gfortran rmdir /s /q build_gfortran
mkdir build_gfortran\mod
set SRC=dependencies\numDeriv-fortran\src\numderiv_kinds.f90 dependencies\numDeriv-fortran\src\numderiv_types.f90 dependencies\numDeriv-fortran\src\numderiv_callbacks.f90 dependencies\numDeriv-fortran\src\numderiv_core.f90 dependencies\numDeriv-fortran\src\numderiv.f90 dependencies\roptim\src\roptim_lbfgsb_core.f90 dependencies\roptim\src\roptim.f90 src\alabama.f90
gfortran -std=f2018 -O0 -g -Wimplicit-interface -Werror=implicit-interface -fcheck=all -J build_gfortran\mod -I build_gfortran\mod %SRC% test\test_alabama.f90 -o build_gfortran\test_alabama.exe
if errorlevel 1 exit /b 1
build_gfortran\test_alabama.exe
if errorlevel 1 exit /b 1
gfortran -std=f2018 -O0 -g -Wimplicit-interface -Werror=implicit-interface -fcheck=all -J build_gfortran\mod -I build_gfortran\mod %SRC% example\constrained_example.f90 -o build_gfortran\constrained_example.exe
if errorlevel 1 exit /b 1
build_gfortran\constrained_example.exe
