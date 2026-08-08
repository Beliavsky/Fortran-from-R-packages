@echo off
setlocal
if exist build-optimized rmdir /s /q build-optimized
mkdir build-optimized\mod
set FLAGS=-std=f2018 -O3 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface
gfortran %FLAGS% -J build-optimized\mod src\mize.f90 test\test_mize.f90 -o build-optimized\test_mize.exe || exit /b 1
build-optimized\test_mize.exe || exit /b 1
gfortran %FLAGS% -J build-optimized\mod src\mize.f90 example\rosenbrock.f90 -o build-optimized\rosenbrock.exe || exit /b 1
build-optimized\rosenbrock.exe || exit /b 1
gfortran %FLAGS% -J build-optimized\mod src\mize.f90 example\stateful.f90 -o build-optimized\stateful.exe || exit /b 1
build-optimized\stateful.exe || exit /b 1
endlocal
