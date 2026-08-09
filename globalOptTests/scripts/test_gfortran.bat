@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%TEMP%\globalOptTests-fortran-build
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\bin"
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0

gfortran %FLAGS% -c "%ROOT%\src\global_opt_tests.f90" -J "%BUILD%\mod" -o "%BUILD%\global_opt_tests.o" || exit /b 1

for %%F in (test_reference_values test_known_optima test_metadata) do (
  gfortran %FLAGS% -I "%BUILD%\mod" "%ROOT%\test\%%F.f90" "%BUILD%\global_opt_tests.o" -o "%BUILD%\bin\%%F.exe" || exit /b 1
  "%BUILD%\bin\%%F.exe" || exit /b 1
)

for %%F in (basic_usage list_benchmarks) do (
  gfortran %FLAGS% -I "%BUILD%\mod" "%ROOT%\example\%%F.f90" "%BUILD%\global_opt_tests.o" -o "%BUILD%\bin\%%F.exe" || exit /b 1
  "%BUILD%\bin\%%F.exe" || exit /b 1
)
endlocal
