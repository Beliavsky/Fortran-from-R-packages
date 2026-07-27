@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build-validation
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod"
mkdir "%BUILD%\bin"
set FLAGS=-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace
set SOURCES=derivmkts_kinds derivmkts_types derivmkts_math derivmkts_rng derivmkts_black_scholes derivmkts_implied derivmkts_bonds derivmkts_asian_analytic derivmkts_barriers derivmkts_perpetual derivmkts_compound derivmkts_jumps derivmkts_binomial derivmkts_simulation derivmkts_asian_mc derivmkts_greeks derivmkts_quincunx derivmkts
set OBJECTS=
for %%S in (%SOURCES%) do (
  gfortran %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" -c "%ROOT%\src\%%S.f90" -o "%BUILD%\%%S.o" || exit /b 1
  set OBJECTS=!OBJECTS! "%BUILD%\%%S.o"
)
ar rcs "%BUILD%\libderivmkts.a" !OBJECTS! || exit /b 1
for %%F in ("%ROOT%\test\*.f90" "%ROOT%\app\*.f90" "%ROOT%\example\*.f90") do (
  set NAME=%%~nF
  gfortran %FLAGS% -I"%BUILD%\mod" "%%F" "%BUILD%\libderivmkts.a" -o "%BUILD%\bin\!NAME!.exe" || exit /b 1
  "%BUILD%\bin\!NAME!.exe" || exit /b 1
)
endlocal
