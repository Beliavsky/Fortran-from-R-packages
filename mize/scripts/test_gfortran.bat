@echo off
setlocal
if exist build-gfortran rmdir /s /q build-gfortran
mkdir build-gfortran\mod
set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
gfortran %FLAGS% -J build-gfortran\mod src\mize.f90 test\test_mize.f90 -o build-gfortran\test_mize.exe || exit /b 1
build-gfortran\test_mize.exe || exit /b 1
endlocal
