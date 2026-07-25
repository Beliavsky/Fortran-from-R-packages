@echo off
setlocal
if "%1"=="" (set MODE=debug) else (set MODE=%1)
if "%FC%"=="" set FC=gfortran
if /I "%MODE%"=="release" (
  set FLAGS=-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace
) else (
  set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -fbacktrace
)
if not exist build\%MODE% mkdir build\%MODE%
cd build\%MODE%
%FC% %FLAGS% -c ..\..\src\fextremes_kinds.f90 ..\..\src\fextremes_rng.f90 ..\..\src\fextremes_stats.f90 ..\..\src\fextremes_optimize.f90 ..\..\src\fextremes_distributions.f90 ..\..\src\fextremes_preprocess.f90 ..\..\src\fextremes_fit.f90 ..\..\src\fextremes_extremal_index.f90 ..\..\src\fextremes_diagnostics.f90 ..\..\src\fextremes_risk.f90 ..\..\src\fextremes_metrics.f90 ..\..\src\fextremes_csv.f90 || exit /b 1
for %%F in (..\..\test\test_*.f90) do (
  %FC% %FLAGS% %%F *.o -o %%~nF.exe || exit /b 1
  %%~nF.exe || exit /b 1
)
echo %MODE% build and tests passed.
