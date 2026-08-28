@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build-gfortran-check
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"
set FLAGS=-std=f2018 -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Wimplicit-interface -Werror=implicit-interface
gfortran %FLAGS% -c "%ROOT%\src\roptim_lbfgsb_core.f90" || exit /b 1
gfortran %FLAGS% -c "%ROOT%\src\roptim.f90" || exit /b 1
gfortran %FLAGS% -c "%ROOT%\test\test_roptim.f90" || exit /b 1
gfortran %FLAGS% -o test_roptim.exe roptim_lbfgsb_core.o roptim.o test_roptim.o || exit /b 1
test_roptim.exe || exit /b 1
