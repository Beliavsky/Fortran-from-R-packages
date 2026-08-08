@echo off
setlocal
if exist build-optimized rmdir /s /q build-optimized
mkdir build-optimized
gfortran -c -std=f2018 -O3 src\lbfgsb3_core.f90 -J build-optimized -o build-optimized\lbfgsb3_core.o || exit /b 1
gfortran -c -std=f2018 -O3 src\lbfgsb3.f90 -I build-optimized -J build-optimized -o build-optimized\lbfgsb3.o || exit /b 1
gfortran -std=f2018 -O3 test\test_lbfgsb3.f90 build-optimized\lbfgsb3.o build-optimized\lbfgsb3_core.o -I build-optimized -J build-optimized -o build-optimized\test_lbfgsb3.exe || exit /b 1
build-optimized\test_lbfgsb3.exe
