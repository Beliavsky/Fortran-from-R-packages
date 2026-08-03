@echo off
setlocal enabledelayedexpansion

set MODE=%1
if "%MODE%"=="" set MODE=check
if /I "%MODE%"=="check" (
  set FLAGS=-std=f2018 -O0 -Wall -Wextra -Werror -fcheck=all -fbacktrace
) else if /I "%MODE%"=="release" (
  set FLAGS=-std=f2018 -O3 -Wall -Wextra -Werror
) else (
  echo Usage: build_all.bat [check^|release]
  exit /b 2
)

cd /d "%~dp0\.."
set BUILD=build\%MODE%
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"

set SOURCES=skellam_kinds skellam_special skellam_distribution skellam_optimization skellam_estimation skellam
set OBJECTS=
for %%S in (%SOURCES%) do (
  gfortran %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" -c "src\%%S.f90" -o "%BUILD%\obj\%%S.o"
  if errorlevel 1 exit /b 1
  set OBJECTS=!OBJECTS! "%BUILD%\obj\%%S.o"
)

for %%F in (test\*.f90 example\*.f90 app\*.f90) do (
  gfortran %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" !OBJECTS! "%%F" -o "%BUILD%\bin\%%~nF.exe"
  if errorlevel 1 exit /b 1
  "%BUILD%\bin\%%~nF.exe"
  if errorlevel 1 exit /b 1
)
