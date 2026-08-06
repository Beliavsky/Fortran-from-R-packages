@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build\gfortran-check
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"

set SOURCES=%ROOT%\src\splines_kinds.f90 %ROOT%\src\splines_linalg.f90 %ROOT%\src\splines_core.f90 %ROOT%\src\splines_basis.f90 %ROOT%\src\splines.f90

gfortran -std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Werror -fcheck=all -fbacktrace -J "%BUILD%" -I "%BUILD%" -o "%BUILD%\test_splines.exe" %SOURCES% "%ROOT%\test\test_splines.f90"
if errorlevel 1 exit /b 1
"%BUILD%\test_splines.exe"
if errorlevel 1 exit /b 1

gfortran -std=f2018 -O3 -Wall -Wextra -Wpedantic -Werror -J "%BUILD%" -I "%BUILD%" -o "%BUILD%\test_splines_optimized.exe" %SOURCES% "%ROOT%\test\test_splines.f90"
if errorlevel 1 exit /b 1
"%BUILD%\test_splines_optimized.exe"
if errorlevel 1 exit /b 1

gfortran -std=f2018 -O2 -Wall -Wextra -Wpedantic -Werror -J "%BUILD%" -I "%BUILD%" -o "%BUILD%\basic_splines.exe" %SOURCES% "%ROOT%\example\basic_splines.f90"
if errorlevel 1 exit /b 1
"%BUILD%\basic_splines.exe"
