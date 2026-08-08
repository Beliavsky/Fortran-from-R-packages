@echo off
setlocal EnableExtensions EnableDelayedExpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build_strict
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
set FFLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
set SRC=calibrar_kinds calibrar_interfaces calibrar_utils calibrar_random calibrar_gradient calibrar_fitness calibrar_splines calibrar_objective calibrar_stopping calibrar_test_functions calibrar_optimization calibrar
set OBJS=
for %%N in (%SRC%) do (
  gfortran %FFLAGS% -c -J "%BUILD%" -I "%BUILD%" "%ROOT%\src\%%N.f90" -o "%BUILD%\%%N.o" || exit /b 1
  set OBJS=!OBJS! "%BUILD%\%%N.o"
)
for %%F in ("%ROOT%\test\*.f90") do (
  set NAME=%%~nF
  gfortran %FFLAGS% -I "%BUILD%" "%%F" !OBJS! -o "%BUILD%\!NAME!.exe" || exit /b 1
  "%BUILD%\!NAME!.exe" || exit /b 1
)
for %%F in ("%ROOT%\example\*.f90") do (
  set NAME=%%~nF
  gfortran %FFLAGS% -I "%BUILD%" "%%F" !OBJS! -o "%BUILD%\!NAME!.exe" || exit /b 1
  "%BUILD%\!NAME!.exe" || exit /b 1
)
endlocal
