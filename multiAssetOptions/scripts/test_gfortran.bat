@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build\gfortran-check
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
set FLAGS=-std=f2018 -Wall -Wextra -Werror -pedantic -fcheck=all -fbacktrace -O0
set SOURCES=mao_kinds mao_status mao_types mao_sparse mao_grid mao_payoff mao_operator mao_pricing multi_asset_options
set OBJECTS=
for %%S in (%SOURCES%) do (
  gfortran %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" -c "%ROOT%\src\%%S.f90" -o "%BUILD%\obj\%%S.o" || exit /b 1
  set OBJECTS=!OBJECTS! "%BUILD%\obj\%%S.o"
)
for %%F in ("%ROOT%\test\*.f90" "%ROOT%\example\*.f90" "%ROOT%\app\*.f90") do (
  set NAME=%%~nF
  gfortran %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" "%%F" !OBJECTS! -o "%BUILD%\bin\!NAME!.exe" || exit /b 1
  "%BUILD%\bin\!NAME!.exe" || exit /b 1
)
endlocal
