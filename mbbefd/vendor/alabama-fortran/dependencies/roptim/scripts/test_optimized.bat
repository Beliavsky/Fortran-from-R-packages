@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build-gfortran-release
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"
set FLAGS=-std=f2018 -O3
gfortran %FLAGS% -c "%ROOT%\src\roptim_lbfgsb_core.f90" || exit /b 1
gfortran %FLAGS% -c "%ROOT%\src\roptim.f90" || exit /b 1
gfortran %FLAGS% -c "%ROOT%\test\test_roptim.f90" || exit /b 1
gfortran %FLAGS% -o test_roptim.exe roptim_lbfgsb_core.o roptim.o test_roptim.o || exit /b 1
test_roptim.exe || exit /b 1
gfortran %FLAGS% -o rosenbrock_methods.exe roptim_lbfgsb_core.o roptim.o "%ROOT%\example\rosenbrock_methods.f90" || exit /b 1
rosenbrock_methods.exe || exit /b 1
gfortran %FLAGS% -o wild_sann.exe roptim_lbfgsb_core.o roptim.o "%ROOT%\example\wild_sann.f90" || exit /b 1
wild_sann.exe || exit /b 1
