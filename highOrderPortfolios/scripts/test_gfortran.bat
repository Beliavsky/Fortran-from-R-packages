@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0\.."

set SOURCES=src\fitheavytail_kinds.f90 src\fitheavytail_status.f90 src\fitheavytail_types.f90 src\fitheavytail_linalg.f90 src\fitheavytail_rng.f90 src\fitheavytail_special.f90 src\fitheavytail_tail.f90 src\fitheavytail_mvst.f90 src\highorder_types.f90 src\highorder_linalg.f90 src\highorder_moments.f90 src\highorder_optimization.f90 src\highorderportfolios.f90
set COMMON=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fbacktrace

call :build debug "%COMMON% -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow" || exit /b 1
call :build release "%COMMON% -O3" || exit /b 1
exit /b 0

:build
set NAME=%~1
set FLAGS=%~2
set DIR=build\%NAME%
if exist "%DIR%" rmdir /s /q "%DIR%"
mkdir "%DIR%\mod"
set OBJECTS=
for %%S in (%SOURCES%) do (
  gfortran %FLAGS% -J "%DIR%\mod" -I "%DIR%\mod" -c %%S -o "%DIR%\%%~nS.o" || exit /b 1
  set OBJECTS=!OBJECTS! "%DIR%\%%~nS.o"
)
for %%S in (test\*.f90) do (
  gfortran %FLAGS% -J "%DIR%\mod" -I "%DIR%\mod" %%S !OBJECTS! -o "%DIR%\%%~nS.exe" || exit /b 1
  "%DIR%\%%~nS.exe" || exit /b 1
)
for %%S in (example\*.f90 app\*.f90) do (
  gfortran %FLAGS% -J "%DIR%\mod" -I "%DIR%\mod" %%S !OBJECTS! -o "%DIR%\%%~nS.exe" || exit /b 1
  "%DIR%\%%~nS.exe" || exit /b 1
)
exit /b 0
