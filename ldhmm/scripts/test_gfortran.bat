REM SPDX-License-Identifier: Artistic-2.0
@echo off
setlocal enabledelayedexpansion

set ROOT=%~dp0..
set BUILD=%ROOT%\build\gfortran-debug
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"

set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -fcheck=all
set SOURCES=ldhmm_kinds ldhmm_status ldhmm_math ldhmm_types ldhmm_distribution ldhmm_parameters ldhmm_modeling ldhmm_simulation ldhmm_optimization ldhmm_series ldhmm
set OBJECTS=

for %%S in (%SOURCES%) do (
    gfortran %FLAGS% -J "%BUILD%" -I "%BUILD%" -c "%ROOT%\src\%%S.f90" -o "%BUILD%\%%S.o"
    if errorlevel 1 exit /b 1
    set OBJECTS=!OBJECTS! "%BUILD%\%%S.o"
)

for %%T in ("%ROOT%\test\*.f90") do (
    echo Running %%~nT
    gfortran %FLAGS% -I "%BUILD%" "%%T" !OBJECTS! -o "%BUILD%\%%~nT.exe"
    if errorlevel 1 exit /b 1
    "%BUILD%\%%~nT.exe"
    if errorlevel 1 exit /b 1
)

echo All GNU Fortran tests passed.
endlocal
