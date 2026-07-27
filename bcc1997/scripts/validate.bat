@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\.validation-build
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -fbacktrace

for %%M in (bcc1997_kinds bcc1997_types bcc1997_quadrature bcc1997_model bcc1997) do (
  gfortran %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" -c "%ROOT%\src\%%M.f90" -o "%BUILD%\obj\%%M.o"
  if errorlevel 1 exit /b 1
)

set OBJS="%BUILD%\obj\bcc1997_kinds.o" "%BUILD%\obj\bcc1997_types.o" "%BUILD%\obj\bcc1997_quadrature.o" "%BUILD%\obj\bcc1997_model.o" "%BUILD%\obj\bcc1997.o"
for %%T in (test_black_scholes_limit test_pricing test_strikes test_validation) do (
  gfortran %FLAGS% -I "%BUILD%\mod" %OBJS% "%ROOT%\test\%%T.f90" -o "%BUILD%\bin\%%T.exe"
  if errorlevel 1 exit /b 1
  "%BUILD%\bin\%%T.exe"
  if errorlevel 1 exit /b 1
)

for %%P in (bcc1997_demo) do (
  gfortran %FLAGS% -I "%BUILD%\mod" %OBJS% "%ROOT%\app\%%P.f90" -o "%BUILD%\bin\%%P.exe"
  if errorlevel 1 exit /b 1
  "%BUILD%\bin\%%P.exe" >nul
  if errorlevel 1 exit /b 1
)
for %%P in (basic_pricing strike_curve) do (
  gfortran %FLAGS% -I "%BUILD%\mod" %OBJS% "%ROOT%\example\%%P.f90" -o "%BUILD%\bin\%%P.exe"
  if errorlevel 1 exit /b 1
  "%BUILD%\bin\%%P.exe" >nul
  if errorlevel 1 exit /b 1
)

rmdir /s /q "%BUILD%"
echo validation: PASS
endlocal
