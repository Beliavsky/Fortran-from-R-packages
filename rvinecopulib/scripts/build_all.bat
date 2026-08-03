@echo off
setlocal enabledelayedexpansion
set MODE=%1
if "%MODE%"=="" set MODE=checked
if /I "%MODE%"=="checked" (
  set BUILD=build_checked
  set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
) else if /I "%MODE%"=="optimized" (
  set BUILD=build_optimized
  set FLAGS=-std=f2018 -O3 -Wall -Wextra -Wimplicit-interface
) else (
  echo usage: scripts\build_all.bat [checked^|optimized]
  exit /b 2
)
cd /d %~dp0\..
if exist %BUILD% rmdir /s /q %BUILD%
mkdir %BUILD%\mod %BUILD%\obj %BUILD%\bin
for %%M in (rvine_kinds rvine_math rvine_bicop rvine_fit rvine_dvine rvine_cvine rvine_tools rvinecopulib) do (
  gfortran %FLAGS% -J %BUILD%\mod -I %BUILD%\mod -c src\%%M.f90 -o %BUILD%\obj\%%M.o || exit /b 1
)
set OBJS=%BUILD%\obj\rvine_kinds.o %BUILD%\obj\rvine_math.o %BUILD%\obj\rvine_bicop.o %BUILD%\obj\rvine_fit.o %BUILD%\obj\rvine_dvine.o %BUILD%\obj\rvine_cvine.o %BUILD%\obj\rvine_tools.o %BUILD%\obj\rvinecopulib.o
for %%F in (test\*.f90 example\*.f90 app\*.f90) do (
  gfortran %FLAGS% -J %BUILD%\mod -I %BUILD%\mod %%F %OBJS% -o %BUILD%\bin\%%~nF.exe || exit /b 1
)
for %%E in (%BUILD%\bin\test_*.exe) do %%E || exit /b 1
for %%E in (%BUILD%\bin\example_*.exe %BUILD%\bin\demo_*.exe) do %%E >nul || exit /b 1
echo %MODE% build: PASS
endlocal
