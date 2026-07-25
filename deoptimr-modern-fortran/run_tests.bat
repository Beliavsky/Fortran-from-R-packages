@echo off
setlocal
if "%FC%"=="" set FC=gfortran
set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -fbacktrace
if not exist build\windows mkdir build\windows
%FC% %FLAGS% -Jbuild\windows -Ibuild\windows -c src\deoptimr_kinds.f90 src\deoptimr_interfaces.f90 src\deoptimr_rng.f90 src\deoptimr_types.f90 src\deoptimr_utils.f90 src\deoptimr_jde.f90 src\deoptimr_ncde.f90 src\deoptimr.f90 test\test_support.f90 || exit /b 1
%FC% %FLAGS% -Jbuild\windows -Ibuild\windows test\test_jde.f90 deoptimr_kinds.o deoptimr_interfaces.o deoptimr_rng.o deoptimr_types.o deoptimr_utils.o deoptimr_jde.o deoptimr_ncde.o deoptimr.o test_support.o -o build\windows\test_jde.exe || exit /b 1
%FC% %FLAGS% -Jbuild\windows -Ibuild\windows test\test_ncde.f90 deoptimr_kinds.o deoptimr_interfaces.o deoptimr_rng.o deoptimr_types.o deoptimr_utils.o deoptimr_jde.o deoptimr_ncde.o deoptimr.o test_support.o -o build\windows\test_ncde.exe || exit /b 1
build\windows\test_jde.exe || exit /b 1
build\windows\test_ncde.exe || exit /b 1
echo Windows tests passed.
