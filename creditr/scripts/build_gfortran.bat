@echo off
setlocal enabledelayedexpansion
set MODE=%1
if "%MODE%"=="" set MODE=debug
if /I "%MODE%"=="debug" (
  set FLAGS=-std=f2018 -Wall -Wextra -Werror -pedantic -fcheck=all -fbacktrace
) else if /I "%MODE%"=="release" (
  set FLAGS=-std=f2018 -O3 -Wall -Wextra -Werror -pedantic
) else (
  echo usage: build_gfortran.bat [debug^|release]
  exit /b 2
)
set BUILD=build\%MODE%
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\obj" "%BUILD%\mod" "%BUILD%\bin"
set SOURCES=src\creditr_kinds.f90 src\creditr_dates.f90 src\creditr_curve.f90 src\creditr_cds.f90 src\creditr_api.f90
set OBJECTS=
for %%F in (%SOURCES%) do (
  gfortran %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" -c %%F -o "%BUILD%\obj\%%~nF.o" || exit /b 1
  set OBJECTS=!OBJECTS! "%BUILD%\obj\%%~nF.o"
)
for %%F in (test\*.f90 app\*.f90 example\*.f90) do (
  gfortran %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" %%F !OBJECTS! -o "%BUILD%\bin\%%~nF.exe" || exit /b 1
  "%BUILD%\bin\%%~nF.exe" || exit /b 1
)
endlocal
