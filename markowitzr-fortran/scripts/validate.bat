@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%TEMP%\markowitzr-fortran-%RANDOM%
mkdir "%BUILD%" || exit /b 1
pushd "%BUILD%" || exit /b 1

set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace -O0 -g
set SOURCES="%ROOT%\src\markowitzr_kinds.f90" "%ROOT%\src\markowitzr_types.f90" "%ROOT%\src\markowitzr_linalg.f90" "%ROOT%\src\markowitzr_moments.f90" "%ROOT%\src\markowitzr_portfolio.f90" "%ROOT%\src\markowitzr.f90"

gfortran %FLAGS% -c %SOURCES% || goto :fail
for %%F in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% *.o "%%F" -o "%%~nF.exe" || goto :fail
  "%%~nF.exe" || goto :fail
)
for %%F in ("%ROOT%\app\*.f90" "%ROOT%\example\*.f90") do (
  gfortran %FLAGS% *.o "%%F" -o "%%~nF.exe" || goto :fail
  "%%~nF.exe" >nul || goto :fail
)
echo validation: PASS
popd
rmdir /s /q "%BUILD%"
exit /b 0

:fail
popd
rmdir /s /q "%BUILD%"
exit /b 1
