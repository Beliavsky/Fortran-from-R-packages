@echo off
setlocal
if not exist build mkdir build
set FFLAGS=-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -fbacktrace -ffree-line-length-none
set SRC=src\fportfolio_kinds.f90 src\fportfolio_types.f90 src\fportfolio_probability.f90 src\fportfolio_linalg.f90 src\fportfolio_statistics.f90 src\fportfolio_risk.f90 src\fportfolio_optimization.f90 src\fportfolio_backtest.f90 src\fportfolio_monitor.f90 src\fportfolio_io.f90
gfortran %FFLAGS% %SRC% test\test_statistics_risk.f90 -llapack -lblas -o build\test_statistics_risk.exe || exit /b 1
gfortran %FFLAGS% %SRC% test\test_optimization.f90 -llapack -lblas -o build\test_optimization.exe || exit /b 1
gfortran %FFLAGS% %SRC% test\test_backtest_monitor.f90 -llapack -lblas -o build\test_backtest_monitor.exe || exit /b 1
build\test_statistics_risk.exe || exit /b 1
build\test_optimization.exe || exit /b 1
build\test_backtest_monitor.exe || exit /b 1
endlocal
