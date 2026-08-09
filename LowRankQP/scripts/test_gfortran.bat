@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build-strict
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0

gfortran %FLAGS% -J "%BUILD%" -I "%BUILD%" -c "%ROOT%\src\lowrankqp_kinds.f90" -o "%BUILD%\kinds.o" || exit /b 1
gfortran %FLAGS% -J "%BUILD%" -I "%BUILD%" -c "%ROOT%\src\lowrankqp_linalg.f90" -o "%BUILD%\linalg.o" || exit /b 1
gfortran %FLAGS% -J "%BUILD%" -I "%BUILD%" -c "%ROOT%\src\lowrankqp.f90" -o "%BUILD%\lowrankqp.o" || exit /b 1
for %%F in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% -J "%BUILD%" -I "%BUILD%" "%BUILD%\kinds.o" "%BUILD%\linalg.o" "%BUILD%\lowrankqp.o" "%%F" -o "%BUILD%\%%~nF.exe" || exit /b 1
  "%BUILD%\%%~nF.exe" || exit /b 1
)
for %%F in ("%ROOT%\example\*.f90") do (
  gfortran %FLAGS% -J "%BUILD%" -I "%BUILD%" "%BUILD%\kinds.o" "%BUILD%\linalg.o" "%BUILD%\lowrankqp.o" "%%F" -o "%BUILD%\%%~nF.exe" || exit /b 1
  "%BUILD%\%%~nF.exe" >nul || exit /b 1
)
echo All strict tests passed.
endlocal
