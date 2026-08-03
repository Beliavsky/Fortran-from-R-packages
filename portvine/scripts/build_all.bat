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
for %%M in (rugarch_kinds rugarch_math rugarch_rng rugarch_gh rugarch_distributions rugarch_optimizer rugarch_linalg rugarch_resampling rugarch_types rugarch_models rugarch_fit rugarch_risk rugarch_arfima rugarch_backtests rugarch_inference rugarch_evaluation rugarch_workflows rugarch_complete rugarch) do (
  gfortran %FLAGS% -J %BUILD%\mod -I %BUILD%\mod -c vendor\rugarch\src\%%M.f90 -o %BUILD%\obj\%%M.o || exit /b 1
)
for %%M in (rvine_kinds rvine_math rvine_bicop rvine_fit rvine_dvine rvine_cvine rvine_tools rvinecopulib) do (
  gfortran %FLAGS% -J %BUILD%\mod -I %BUILD%\mod -c vendor\rvinecopulib\src\%%M.f90 -o %BUILD%\obj\%%M.o || exit /b 1
)
for %%M in (portvine_kinds portvine_types portvine_stats portvine_ordering portvine_marginals portvine_conditional portvine_dependence portvine_workflow portvine) do (
  gfortran %FLAGS% -J %BUILD%\mod -I %BUILD%\mod -c src\%%M.f90 -o %BUILD%\obj\%%M.o || exit /b 1
)
set OBJS=
for %%O in (%BUILD%\obj\*.o) do set OBJS=!OBJS! %%O
for %%F in (test\*.f90 example\*.f90 app\*.f90) do (
  gfortran %FLAGS% -J %BUILD%\mod -I %BUILD%\mod %%F !OBJS! -o %BUILD%\bin\%%~nF.exe || exit /b 1
)
for %%E in (%BUILD%\bin\test_*.exe) do %%E || exit /b 1
for %%E in (%BUILD%\bin\example_*.exe %BUILD%\bin\demo_*.exe) do %%E >nul || exit /b 1
echo %MODE% build: PASS
endlocal
