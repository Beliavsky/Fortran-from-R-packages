@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build-windows
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
set FLAGS=-std=f2018 -O2 -Wall -Wextra -Werror -J "%BUILD%\mod" -I "%BUILD%\mod"
set SOURCES=kind_mod cla_types cla_core cla_queries cla_garch cla
set OBJECTS=
for %%S in (%SOURCES%) do (
  gfortran %FLAGS% -c "%ROOT%\src\%%S.f90" -o "%BUILD%\obj\%%S.o" || exit /b 1
  set OBJECTS=!OBJECTS! "%BUILD%\obj\%%S.o"
)
for %%F in ("%ROOT%\test\*.f90" "%ROOT%\app\*.f90" "%ROOT%\example\*.f90") do (
  gfortran %FLAGS% "%%F" !OBJECTS! -llapack -lblas -o "%BUILD%\bin\%%~nF.exe" || exit /b 1
  "%BUILD%\bin\%%~nF.exe" || exit /b 1
)
echo CLA Windows build and tests passed
