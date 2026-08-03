@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0\.."
set BUILD=build\optimized
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\bin"
set FLAGS=-O3 -std=f2018 -Wall -Wextra -Werror -Wimplicit-interface
set SOURCES=vendor\rugarch-modern-fortran\src\rugarch_kinds.f90 vendor\rugarch-modern-fortran\src\rugarch_math.f90 vendor\rugarch-modern-fortran\src\rugarch_rng.f90 vendor\rugarch-modern-fortran\src\rugarch_gh.f90 vendor\rugarch-modern-fortran\src\rugarch_distributions.f90 vendor\rugarch-modern-fortran\src\rugarch_optimizer.f90 vendor\rugarch-modern-fortran\src\rugarch_linalg.f90 vendor\rugarch-modern-fortran\src\rugarch_resampling.f90 vendor\rugarch-modern-fortran\src\rugarch_types.f90 vendor\rugarch-modern-fortran\src\rugarch_models.f90 vendor\rugarch-modern-fortran\src\rugarch_fit.f90 vendor\rugarch-modern-fortran\src\rugarch_risk.f90 vendor\rugarch-modern-fortran\src\rugarch_arfima.f90 vendor\rugarch-modern-fortran\src\rugarch_backtests.f90 vendor\rugarch-modern-fortran\src\rugarch_inference.f90 vendor\rugarch-modern-fortran\src\rugarch_evaluation.f90 vendor\rugarch-modern-fortran\src\rugarch_workflows.f90 vendor\rugarch-modern-fortran\src\rugarch_complete.f90 vendor\rugarch-modern-fortran\src\rugarch.f90 src\quarks_kinds.f90 src\quarks_types.f90 src\quarks_rng.f90 src\quarks_stats.f90 src\quarks_smoothing.f90 src\quarks_risk.f90 src\quarks_backtests.f90 src\quarks_portfolio.f90 src\quarks_rollcast.f90 src\quarks.f90
for %%S in (%SOURCES%) do gfortran %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" -c %%S -o "%BUILD%\%%~nS.o" || exit /b 1
ar rcs "%BUILD%\libquarks.a" "%BUILD%\*.o" || exit /b 1
for %%T in (test\*.f90) do (
  gfortran %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" %%T "%BUILD%\libquarks.a" -o "%BUILD%\bin\%%~nT.exe" || exit /b 1
  "%BUILD%\bin\%%~nT.exe" || exit /b 1
)
for %%T in (example\*.f90 app\*.f90) do (
  gfortran %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" %%T "%BUILD%\libquarks.a" -o "%BUILD%\bin\%%~nT.exe" || exit /b 1
  "%BUILD%\bin\%%~nT.exe" || exit /b 1
)
