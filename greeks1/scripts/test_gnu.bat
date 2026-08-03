@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build-gfortran
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow
set SOURCES=greeks_kinds greeks_types greeks_math greeks_payoffs greeks_integrals greeks_black_scholes greeks_binomial greeks_monte_carlo greeks_api greeks
set OBJECTS=
for %%S in (%SOURCES%) do (
  gfortran %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" -c "%ROOT%\src\%%S.f90" -o "%BUILD%\obj\%%S.o" || exit /b 1
  set OBJECTS=!OBJECTS! "%BUILD%\obj\%%S.o"
)
for %%T in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" "%%T" !OBJECTS! -o "%BUILD%\bin\%%~nT.exe" || exit /b 1
  "%BUILD%\bin\%%~nT.exe" || exit /b 1
)
endlocal
