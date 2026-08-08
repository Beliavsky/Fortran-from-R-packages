@echo off
setlocal
if exist build-debug rmdir /s /q build-debug
mkdir build-debug
gfortran -std=f2018 -Wall -Wextra -Wpedantic -Wno-compare-reals -Wimplicit-interface -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -J build-debug -I build-debug src\n1qn1.f90 test\test_n1qn1.f90 -o build-debug\test_n1qn1.exe
if errorlevel 1 exit /b 1
build-debug\test_n1qn1.exe
if errorlevel 1 exit /b 1
gfortran -std=f2018 -Wall -Wextra -Wpedantic -Wno-compare-reals -Wimplicit-interface -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -J build-debug -I build-debug src\n1qn1.f90 example\banana.f90 -o build-debug\banana.exe
if errorlevel 1 exit /b 1
build-debug\banana.exe
