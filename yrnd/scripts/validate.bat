@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build-validation-windows
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"
set TVM=%ROOT%\dependencies\tvm-fortran\src
set FLAGS=-std=f2018 -O2 -Wall -Wextra -fcheck=all -fbacktrace -J. -I.
gfortran %FLAGS% -c "%TVM%\tvm_kinds.f90" "%TVM%\tvm_root.f90" "%TVM%\tvm_interpolation.f90" "%TVM%\tvm_cashflows.f90" "%TVM%\tvm_curves.f90" "%TVM%\tvm.f90" "%ROOT%\src\yrnd_kinds.f90" "%ROOT%\src\yrnd_dates.f90" "%ROOT%\src\yrnd_stats.f90" "%ROOT%\src\yrnd_optimize.f90" "%ROOT%\src\yrnd_mixture.f90" "%ROOT%\src\yrnd_bonds.f90" "%ROOT%\src\yrnd_transforms.f90" "%ROOT%\src\yrnd_api.f90" "%ROOT%\src\yrnd.f90"
if errorlevel 1 exit /b 1
gfortran %FLAGS% -c "%ROOT%\test\test_yrnd.f90"
gfortran %FLAGS% -o test_yrnd.exe *.o
if errorlevel 1 exit /b 1
test_yrnd.exe
if errorlevel 1 exit /b 1
del test_yrnd.o test_yrnd.exe
gfortran %FLAGS% -c "%ROOT%\test\test_api.f90"
gfortran %FLAGS% -o test_api.exe *.o
if errorlevel 1 exit /b 1
test_api.exe
