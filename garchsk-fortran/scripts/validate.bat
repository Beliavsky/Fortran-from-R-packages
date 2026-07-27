@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build-validation
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"
set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Werror -Wconversion-extra -Wimplicit-interface -fcheck=all -fbacktrace -J . -I .
for %%F in (garchsk_kinds garchsk_types garchsk_stats garchsk_linalg garchsk_models garchsk_estimation garchsk) do (
  gfortran %FLAGS% -c "%ROOT%\src\%%F.f90" || exit /b 1
)
gfortran %FLAGS% -c "%ROOT%\test\test_support.f90" || exit /b 1
set LIB=garchsk_kinds.o garchsk_types.o garchsk_stats.o garchsk_linalg.o garchsk_models.o garchsk_estimation.o garchsk.o
for %%T in (test_constraints test_estimation test_likelihood_forecast test_moments) do (
  gfortran %FLAGS% "%ROOT%\test\%%T.f90" %LIB% test_support.o -o %%T.exe || exit /b 1
  %%T.exe || exit /b 1
)
for %%P in (garchsk_demo basic_garchsk estimate_and_forecast) do (
  if "%%P"=="garchsk_demo" set SRC=%ROOT%\app\%%P.f90
  if not "%%P"=="garchsk_demo" set SRC=%ROOT%\example\%%P.f90
  gfortran %FLAGS% "!SRC!" %LIB% -o %%P.exe || exit /b 1
  %%P.exe >nul || exit /b 1
)
echo validation: PASS
endlocal
