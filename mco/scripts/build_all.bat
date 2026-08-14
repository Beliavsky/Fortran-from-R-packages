@echo off
setlocal enabledelayedexpansion
set MODE=%1
if "%MODE%"=="" set MODE=checked
if "%MODE%"=="checked" (
  set FLAGS=-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all -fbacktrace -O0
) else if "%MODE%"=="optimized" (
  set FLAGS=-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -O3
) else (
  echo usage: build_all.bat checked^|optimized
  exit /b 2
)
set BUILD=build\%MODE%
if exist %BUILD% rmdir /s /q %BUILD%
mkdir %BUILD%\mod %BUILD%\obj %BUILD%\bin
set SOURCES=mco_kinds mco_random mco_pareto mco_quality mco_nsga2 mco_test_functions mco
set OBJECTS=
for %%S in (%SOURCES%) do (
  gfortran %FLAGS% -J%BUILD%\mod -I%BUILD%\mod -c src\%%S.f90 -o %BUILD%\obj\%%S.o || exit /b 1
  set OBJECTS=!OBJECTS! %BUILD%\obj\%%S.o
)
for %%F in (test\*.f90 example\*.f90 app\*.f90) do (
  gfortran %FLAGS% -I%BUILD%\mod %%F !OBJECTS! -o %BUILD%\bin\%%~nF.exe || exit /b 1
  %BUILD%\bin\%%~nF.exe || exit /b 1
)
