@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build-gfortran
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"

set FLAGS=-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
set SOURCES="%ROOT%\src\dfoptim_kinds.f90" "%ROOT%\src\dfoptim_interfaces.f90" "%ROOT%\src\dfoptim_rng.f90" "%ROOT%\src\dfoptim_utils.f90" "%ROOT%\src\dfoptim_hooke_jeeves.f90" "%ROOT%\src\dfoptim_nelder_mead.f90" "%ROOT%\src\dfoptim_mads.f90" "%ROOT%\src\dfoptim.f90"

gfortran %FLAGS% -J "%BUILD%" -I "%BUILD%" %SOURCES% "%ROOT%\test\test_dfoptim.f90" -o "%BUILD%\test_dfoptim.exe" || exit /b 1
"%BUILD%\test_dfoptim.exe" || exit /b 1

gfortran %FLAGS% -J "%BUILD%" -I "%BUILD%" %SOURCES% "%ROOT%\example\dfoptim_example.f90" -o "%BUILD%\dfoptim_example.exe" || exit /b 1
"%BUILD%\dfoptim_example.exe" || exit /b 1
endlocal
