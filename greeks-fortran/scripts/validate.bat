@echo off
setlocal enabledelayedexpansion
if "%FC%"=="" set FC=gfortran
set ROOT=%~dp0..
set SOURCES=%ROOT%\src\greeks_kinds.f90 %ROOT%\src\greeks_types.f90 %ROOT%\src\greeks_math.f90 %ROOT%\src\greeks_rng.f90 %ROOT%\src\greeks_european.f90 %ROOT%\src\greeks_geometric_asian.f90 %ROOT%\src\greeks_american.f90 %ROOT%\src\greeks_paths.f90 %ROOT%\src\greeks_malliavin.f90 %ROOT%\src\greeks_implied_volatility.f90 %ROOT%\src\greeks.f90

call :build debug "-std=f2018 -O0 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace"
if errorlevel 1 exit /b 1
call :build release "-std=f2018 -O2 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror"
if errorlevel 1 exit /b 1
echo validation: PASS
exit /b 0

:build
set MODE=%~1
set FLAGS=%~2
set BUILD=%ROOT%\build\%MODE%
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
pushd "%BUILD%"
%FC% %FLAGS% -J. -I. -c %SOURCES%
if errorlevel 1 goto :failed
for %%F in ("%ROOT%\test\*.f90") do (
  set NAME=%%~nF
  %FC% %FLAGS% -J. -I. *.o "%%F" -o !NAME!.exe
  if errorlevel 1 goto :failed
  !NAME!.exe
  if errorlevel 1 goto :failed
)
for %%F in ("%ROOT%\app\*.f90" "%ROOT%\example\*.f90") do (
  set NAME=%%~nF
  %FC% %FLAGS% -J. -I. *.o "%%F" -o !NAME!.exe
  if errorlevel 1 goto :failed
  !NAME!.exe >nul
  if errorlevel 1 goto :failed
)
popd
exit /b 0
:failed
popd
exit /b 1
