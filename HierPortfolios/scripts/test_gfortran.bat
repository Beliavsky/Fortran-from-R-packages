@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build-gfortran
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow

gfortran %FLAGS% -J . -I . -c ^
  "%ROOT%\src\hierportfolios_kinds.f90" ^
  "%ROOT%\src\hierportfolios_types.f90" ^
  "%ROOT%\src\hierportfolios_hierarchy.f90" ^
  "%ROOT%\src\hierportfolios_core.f90" ^
  "%ROOT%\src\hierportfolios.f90" || exit /b 1

for %%F in ("%ROOT%\test\*.f90" "%ROOT%\example\*.f90" "%ROOT%\app\*.f90") do (
  gfortran %FLAGS% -J . -I . *.o "%%~fF" -o "%%~nF.exe" || exit /b 1
  "%%~nF.exe" || exit /b 1
)
endlocal
