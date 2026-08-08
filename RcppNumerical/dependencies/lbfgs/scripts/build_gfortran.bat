@echo off
setlocal
set BUILD=build-gfortran
if not "%~1"=="" set BUILD=%~1
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"

set FLAGS=-std=f2018 -O2 -Wall -Wextra -Wpedantic
set MODS=-J "%BUILD%" -I "%BUILD%"
set SRC=src\lbfgs_kinds.f90 src\lbfgs_status.f90 src\lbfgs_solver.f90 src\lbfgs.f90

gfortran %FLAGS% %MODS% %SRC% test\test_lbfgs.f90 -o "%BUILD%\test_lbfgs.exe" || exit /b 1
"%BUILD%\test_lbfgs.exe" || exit /b 1

gfortran %FLAGS% %MODS% %SRC% example\rosenbrock.f90 -o "%BUILD%\rosenbrock.exe" || exit /b 1
"%BUILD%\rosenbrock.exe" || exit /b 1

gfortran %FLAGS% %MODS% %SRC% example\owlqn_soft_threshold.f90 -o "%BUILD%\owlqn_soft_threshold.exe" || exit /b 1
"%BUILD%\owlqn_soft_threshold.exe" || exit /b 1
endlocal
