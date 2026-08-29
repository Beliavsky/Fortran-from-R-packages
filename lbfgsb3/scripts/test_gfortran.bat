@echo off
setlocal
if exist build-gfortran rmdir /s /q build-gfortran
mkdir build-gfortran
gfortran -c -std=f2018 -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow src\lbfgsb3_core.f90 -J build-gfortran -o build-gfortran\lbfgsb3_core.o || exit /b 1
gfortran -c -std=f2018 -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow src\lbfgsb3.f90 -I build-gfortran -J build-gfortran -o build-gfortran\lbfgsb3.o || exit /b 1
gfortran -std=f2018 -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow test\test_lbfgsb3.f90 build-gfortran\lbfgsb3.o build-gfortran\lbfgsb3_core.o -I build-gfortran -J build-gfortran -o build-gfortran\test_lbfgsb3.exe || exit /b 1
build-gfortran\test_lbfgsb3.exe
