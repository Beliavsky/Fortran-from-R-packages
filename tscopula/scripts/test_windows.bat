@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build\windows-check
set MOD=%BUILD%\mod
if not exist "%BUILD%" mkdir "%BUILD%"
if not exist "%MOD%" mkdir "%MOD%"
set FLAGS=-std=f2018 -ffree-line-length-none -Wall -Wextra -Werror -Wno-maybe-uninitialized -O0 -g -fcheck=all -fbacktrace -J"%MOD%" -I"%MOD%"
set SOURCES=tscopula_kinds tscopula_status tscopula_math tscopula_margins tscopula_vtransforms tscopula_paircopula tscopula_timeseries tscopula_dvine tscopula_models tscopula_compat tscopula
set OBJECTS=
for %%S in (%SOURCES%) do (
  gfortran %FLAGS% -c "%ROOT%\src\%%S.f90" -o "%BUILD%\%%S.o" || exit /b 1
  set OBJECTS=!OBJECTS! "%BUILD%\%%S.o"
)
gfortran %FLAGS% -c "%ROOT%\test\test_utils.f90" -o "%BUILD%\test_utils.o" || exit /b 1
set OBJECTS=%OBJECTS% "%BUILD%\test_utils.o"
for %%T in (test_margins_vtrans test_paircopula test_timeseries test_dvine test_models test_fitting) do (
  gfortran %FLAGS% "%ROOT%\test\%%T.f90" %OBJECTS% -o "%BUILD%\%%T.exe" || exit /b 1
  "%BUILD%\%%T.exe" || exit /b 1
)
endlocal
