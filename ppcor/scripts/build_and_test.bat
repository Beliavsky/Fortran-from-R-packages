@echo off
setlocal enabledelayedexpansion

if "%FC%"=="" set FC=gfortran
set MODE=%1
if "%MODE%"=="" set MODE=checked
set ROOT=%~dp0..
set BUILD=%ROOT%\build-manual-%MODE%

if /I "%MODE%"=="checked" (
  set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Wconversion-extra -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow
) else if /I "%MODE%"=="optimized" (
  set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Wconversion-extra -O3
) else (
  echo Usage: scripts\build_and_test.bat [checked^|optimized]
  exit /b 2
)

if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin" || exit /b 1

call :compile src\ppcor_kinds.f90 ppcor_kinds.o || exit /b 1
call :compile src\ppcor_special.f90 ppcor_special.o || exit /b 1
call :compile src\ppcor_linalg.f90 ppcor_linalg.o || exit /b 1
call :compile src\ppcor_stats.f90 ppcor_stats.o || exit /b 1
call :compile src\ppcor.f90 ppcor.o || exit /b 1

set OBJS="%BUILD%\obj\ppcor_kinds.o" "%BUILD%\obj\ppcor_special.o" "%BUILD%\obj\ppcor_linalg.o" "%BUILD%\obj\ppcor_stats.o" "%BUILD%\obj\ppcor.o"

for %%F in ("%ROOT%\test\*.f90") do (
  %FC% %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" %OBJS% "%%F" -o "%BUILD%\bin\%%~nF.exe" || exit /b 1
  "%BUILD%\bin\%%~nF.exe" || exit /b 1
)

for %%F in ("%ROOT%\example\*.f90" "%ROOT%\app\*.f90") do (
  %FC% %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" %OBJS% "%%F" -o "%BUILD%\bin\%%~nF.exe" || exit /b 1
  "%BUILD%\bin\%%~nF.exe" >nul || exit /b 1
)

echo All %MODE% tests and examples passed.
exit /b 0

:compile
%FC% %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" -c "%ROOT%\%~1" -o "%BUILD%\obj\%~2"
exit /b %errorlevel%
