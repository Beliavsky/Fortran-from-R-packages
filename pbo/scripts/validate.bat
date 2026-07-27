@echo off
setlocal
cd /d %~dp0\..
if exist build-validation rmdir /s /q build-validation
mkdir build-validation
set FLAGS=-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace -Jbuild-validation -Ibuild-validation

gfortran %FLAGS% -c src\pbo_kinds.f90 -o build-validation\pbo_kinds.o || exit /b 1
gfortran %FLAGS% -c src\pbo_types.f90 -o build-validation\pbo_types.o || exit /b 1
gfortran %FLAGS% -c src\pbo_combinations.f90 -o build-validation\pbo_combinations.o || exit /b 1
gfortran %FLAGS% -c src\pbo_stats.f90 -o build-validation\pbo_stats.o || exit /b 1
gfortran %FLAGS% -c src\pbo_metrics.f90 -o build-validation\pbo_metrics.o || exit /b 1
gfortran %FLAGS% -c src\pbo_core.f90 -o build-validation\pbo_core.o || exit /b 1
gfortran %FLAGS% -c src\pbo_analysis.f90 -o build-validation\pbo_analysis.o || exit /b 1
gfortran %FLAGS% -c src\pbo.f90 -o build-validation\pbo.o || exit /b 1
set OBJS=build-validation\pbo_kinds.o build-validation\pbo_types.o build-validation\pbo_combinations.o build-validation\pbo_stats.o build-validation\pbo_metrics.o build-validation\pbo_core.o build-validation\pbo_analysis.o build-validation\pbo.o

for %%F in (test\*.f90 app\*.f90 example\*.f90) do (
  gfortran %FLAGS% %%F %OBJS% -o build-validation\%%~nF.exe || exit /b 1
  build-validation\%%~nF.exe || exit /b 1
)
echo validation: PASS
endlocal
