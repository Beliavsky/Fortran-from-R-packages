@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build_strict
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all

gfortran %FLAGS% -c "%ROOT%\src\neighbours_kinds.f90" "%ROOT%\src\neighbours_rng.f90" "%ROOT%\src\neighbours.f90" || exit /b 1
for %%F in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% "%%F" neighbours_kinds.o neighbours_rng.o neighbours.o -o "%%~nF.exe" || exit /b 1
  "%%~nF.exe" || exit /b 1
)
for %%F in ("%ROOT%\example\*.f90") do (
  gfortran %FLAGS% "%%F" neighbours_kinds.o neighbours_rng.o neighbours.o -o "%%~nF.exe" || exit /b 1
  "%%~nF.exe" || exit /b 1
)
endlocal
