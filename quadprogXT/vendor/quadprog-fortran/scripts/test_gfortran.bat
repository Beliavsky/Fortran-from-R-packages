@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build-gfortran-windows
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Werror -Wno-compare-reals -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow

gfortran %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" -c "%ROOT%\src\quadprog_kinds.f90" -o "%BUILD%\obj\quadprog_kinds.o" || exit /b 1
gfortran %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" -c "%ROOT%\src\quadprog_core.f90" -o "%BUILD%\obj\quadprog_core.o" || exit /b 1
gfortran %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" -c "%ROOT%\src\quadprog.f90" -o "%BUILD%\obj\quadprog.o" || exit /b 1
gfortran %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" -c "%ROOT%\src\test_support.f90" -o "%BUILD%\obj\test_support.o" || exit /b 1

for %%F in ("%ROOT%\test\*.f90") do (
  set NAME=%%~nF
  gfortran %FLAGS% -I"%BUILD%\mod" "%%F" "%BUILD%\obj\quadprog_kinds.o" "%BUILD%\obj\quadprog_core.o" "%BUILD%\obj\quadprog.o" "%BUILD%\obj\test_support.o" -o "%BUILD%\bin\!NAME!.exe" || exit /b 1
  "%BUILD%\bin\!NAME!.exe" || exit /b 1
)

echo GNU Fortran Windows validation: PASS
endlocal
