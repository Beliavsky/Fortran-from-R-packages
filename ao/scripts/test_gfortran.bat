@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build-gfortran
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
gfortran %FLAGS% -J. -I. -c "%ROOT%\src\ao_kinds.f90" "%ROOT%\src\ao_types.f90" "%ROOT%\src\ao_random.f90" "%ROOT%\src\ao_history.f90" "%ROOT%\src\ao_base_optimizer.f90" "%ROOT%\src\ao.f90" || exit /b 1
for %%F in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% -I. *.o "%%F" -o test.exe || exit /b 1
  test.exe || exit /b 1
)
for %%F in ("%ROOT%\example\*.f90") do (
  gfortran %FLAGS% -I. *.o "%%F" -o example.exe || exit /b 1
  example.exe || exit /b 1
)
endlocal
