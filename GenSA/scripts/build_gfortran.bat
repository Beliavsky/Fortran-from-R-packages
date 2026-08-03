@echo off
setlocal enabledelayedexpansion
set MODE=%1
if "%MODE%"=="" set MODE=checked
set ROOT=%~dp0..
set BUILD=%ROOT%\build-%MODE%
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"

set COMMON=-std=f2018 -Wall -Wextra -Wpedantic -Wconversion-extra
if /I "%MODE%"=="checked" (
  set FLAGS=-O0 -g -fcheck=all -fbacktrace
) else if /I "%MODE%"=="optimized" (
  set FLAGS=-O3
) else (
  echo Usage: build_gfortran.bat [checked^|optimized]
  exit /b 2
)

set OBJECTS=
for %%S in (gensa_kinds gensa_rng gensa_types gensa_local gensa) do (
  gfortran %COMMON% %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" -c "%ROOT%\src\%%S.f90" -o "%BUILD%\obj\%%S.o" || exit /b 1
  set OBJECTS=!OBJECTS! "%BUILD%\obj\%%S.o"
)

for %%F in ("%ROOT%\test\*.f90") do (
  gfortran %COMMON% %FLAGS% -I"%BUILD%\mod" "%%F" !OBJECTS! -o "%BUILD%\bin\%%~nF.exe" || exit /b 1
  "%BUILD%\bin\%%~nF.exe" || exit /b 1
)

for %%F in ("%ROOT%\example\*.f90" "%ROOT%\app\*.f90") do (
  gfortran %COMMON% %FLAGS% -I"%BUILD%\mod" "%%F" !OBJECTS! -o "%BUILD%\bin\%%~nF.exe" || exit /b 1
  "%BUILD%\bin\%%~nF.exe" >nul || exit /b 1
)

echo All %MODE% tests and examples passed.
endlocal
