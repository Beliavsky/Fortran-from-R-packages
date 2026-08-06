@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build\tests
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"
set SRC=%ROOT%\dependencies\splines\src\splines_kinds.f90 %ROOT%\dependencies\splines\src\splines_linalg.f90 %ROOT%\dependencies\splines\src\splines_core.f90 %ROOT%\dependencies\splines\src\splines_basis.f90 %ROOT%\dependencies\splines\src\splines.f90 %ROOT%\src\mgcv_kinds.f90 %ROOT%\src\mgcv_linalg.f90 %ROOT%\src\mgcv_utils.f90 %ROOT%\src\mgcv_distributions.f90 %ROOT%\src\mgcv_smooths.f90 %ROOT%\src\mgcv_families.f90 %ROOT%\src\mgcv_fit.f90 %ROOT%\src\mgcv_constraints.f90 %ROOT%\src\mgcv_discrete.f90 %ROOT%\src\mgcv_matrix.f90 %ROOT%\src\mgcv_simulation.f90 %ROOT%\src\mgcv.f90
gfortran -std=f2018 -Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace -ffree-line-length-none -O0 -g %SRC% "%ROOT%\test\test_mgcv.f90" -o test_mgcv.exe
if errorlevel 1 exit /b 1
test_mgcv.exe
