@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0\.."
if exist build-gfortran-optimized rmdir /s /q build-gfortran-optimized
mkdir build-gfortran-optimized\mod build-gfortran-optimized\obj build-gfortran-optimized\bin
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -Wimplicit-interface -O3
set SOURCES=vendor\rumidas-fortran\vendor\maxLik-fortran\src\maxlik_kinds.f90 vendor\rumidas-fortran\vendor\maxLik-fortran\src\maxlik_status.f90 vendor\rumidas-fortran\vendor\maxLik-fortran\src\maxlik_types.f90 vendor\rumidas-fortran\vendor\maxLik-fortran\src\maxlik_linalg.f90 vendor\rumidas-fortran\vendor\maxLik-fortran\src\maxlik_random.f90 vendor\rumidas-fortran\vendor\maxLik-fortran\src\maxlik_evaluation.f90 vendor\rumidas-fortran\vendor\maxLik-fortran\src\maxlik_solvers.f90 vendor\rumidas-fortran\vendor\maxLik-fortran\src\maxlik_inference.f90 vendor\rumidas-fortran\vendor\maxLik-fortran\src\maxlik_utilities.f90 vendor\rumidas-fortran\vendor\maxLik-fortran\src\maxlik_api.f90 vendor\rumidas-fortran\vendor\maxLik-fortran\src\maxlik.f90 vendor\rumidas-fortran\src\rumidas_kinds.f90 vendor\rumidas-fortran\src\rumidas_status.f90 vendor\rumidas-fortran\src\rumidas_types.f90 vendor\rumidas-fortran\src\rumidas_weights.f90 vendor\rumidas-fortran\src\rumidas_statistics.f90 vendor\rumidas-fortran\src\rumidas_garch_midas.f90 vendor\rumidas-fortran\src\rumidas_mem.f90 vendor\rumidas-fortran\src\rumidas_fit.f90 vendor\rumidas-fortran\src\rumidas_forecast.f90 vendor\rumidas-fortran\src\rumidas.f90 vendor\rugarch-modern-fortran\src\rugarch_kinds.f90 vendor\rugarch-modern-fortran\src\rugarch_math.f90 vendor\rugarch-modern-fortran\src\rugarch_rng.f90 vendor\rugarch-modern-fortran\src\rugarch_gh.f90 vendor\rugarch-modern-fortran\src\rugarch_distributions.f90 vendor\rugarch-modern-fortran\src\rugarch_optimizer.f90 vendor\rugarch-modern-fortran\src\rugarch_linalg.f90 vendor\rugarch-modern-fortran\src\rugarch_resampling.f90 vendor\rugarch-modern-fortran\src\rugarch_types.f90 vendor\rugarch-modern-fortran\src\rugarch_models.f90 vendor\rugarch-modern-fortran\src\rugarch_fit.f90 vendor\rugarch-modern-fortran\src\rugarch_risk.f90 vendor\rugarch-modern-fortran\src\rugarch_arfima.f90 vendor\rugarch-modern-fortran\src\rugarch_backtests.f90 vendor\rugarch-modern-fortran\src\rugarch_inference.f90 vendor\rugarch-modern-fortran\src\rugarch_evaluation.f90 vendor\rugarch-modern-fortran\src\rugarch_workflows.f90 vendor\rugarch-modern-fortran\src\rugarch_complete.f90 vendor\rugarch-modern-fortran\src\rugarch.f90 src\pwev_kinds.f90 src\pwev_status.f90 src\pwev_types.f90 src\pwev_metrics.f90 src\pwev_pso.f90 src\pwev_models.f90 src\pwev_core.f90 src\pwev.f90
for %%F in (%SOURCES%) do (
  gfortran %FLAGS% -Jbuild-gfortran-optimized\mod -Ibuild-gfortran-optimized\mod -c %%F -o build-gfortran-optimized\obj\%%~nF.o || exit /b 1
)
ar rcs build-gfortran-optimized\libpwev.a build-gfortran-optimized\obj\*.o || exit /b 1
for %%F in (test\*.f90) do (
  gfortran %FLAGS% -Ibuild-gfortran-optimized\mod %%F build-gfortran-optimized\libpwev.a -o build-gfortran-optimized\bin\%%~nF.exe || exit /b 1
  build-gfortran-optimized\bin\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90 app\*.f90) do (
  gfortran %FLAGS% -Ibuild-gfortran-optimized\mod %%F build-gfortran-optimized\libpwev.a -o build-gfortran-optimized\bin\%%~nF.exe || exit /b 1
  build-gfortran-optimized\bin\%%~nF.exe >nul || exit /b 1
)
echo optimized GNU Fortran validation: PASS
