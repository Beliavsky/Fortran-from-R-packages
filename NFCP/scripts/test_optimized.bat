@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build-gfortran-opt
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"
set FLAGS=-std=f2018 -O3 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface
gfortran %FLAGS% -J . -I . -c ^
 "%ROOT%\dependencies\lbfgsb3\src\lbfgsb3_core.f90" ^
 "%ROOT%\dependencies\lbfgsb3\src\lbfgsb3.f90" ^
 "%ROOT%\src\nfcp_types.f90" "%ROOT%\src\nfcp_math.f90" ^
 "%ROOT%\src\nfcp_parameters.f90" "%ROOT%\src\nfcp_kalman.f90" ^
 "%ROOT%\src\nfcp_forecast.f90" "%ROOT%\src\nfcp_simulation.f90" ^
 "%ROOT%\src\nfcp_analysis.f90" "%ROOT%\src\nfcp_stitch.f90" ^
 "%ROOT%\src\nfcp_options.f90" "%ROOT%\src\nfcp_mle.f90" "%ROOT%\src\nfcp.f90"
if errorlevel 1 exit /b 1
gfortran %FLAGS% -I . "%ROOT%\test\test_nfcp.f90" *.o -o test_nfcp.exe
if errorlevel 1 exit /b 1
test_nfcp.exe
if errorlevel 1 exit /b 1
gfortran %FLAGS% -I . "%ROOT%\example\two_factor_oil.f90" *.o -o two_factor_oil.exe
if errorlevel 1 exit /b 1
two_factor_oil.exe
endlocal
