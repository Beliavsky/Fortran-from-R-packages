@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build\example
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"
set SRC=%ROOT%\dependencies\splines\src\splines_kinds.f90 %ROOT%\dependencies\splines\src\splines_linalg.f90 %ROOT%\dependencies\splines\src\splines_core.f90 %ROOT%\dependencies\splines\src\splines_basis.f90 %ROOT%\dependencies\splines\src\splines.f90 %ROOT%\src\mgcv_kinds.f90 %ROOT%\src\mgcv_linalg.f90 %ROOT%\src\mgcv_utils.f90 %ROOT%\src\mgcv_distributions.f90 %ROOT%\src\mgcv_smooths.f90 %ROOT%\src\mgcv_families.f90 %ROOT%\src\mgcv_fit.f90 %ROOT%\src\mgcv_constraints.f90 %ROOT%\src\mgcv_discrete.f90 %ROOT%\src\mgcv_matrix.f90 %ROOT%\src\mgcv_simulation.f90 %ROOT%\src\mgcv.f90
gfortran -std=f2018 -O2 -ffree-line-length-none %SRC% "%ROOT%\example\basic_gam.f90" -o basic_gam.exe
if errorlevel 1 exit /b 1
basic_gam.exe
