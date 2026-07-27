@echo off
setlocal enabledelayedexpansion
set MODE=%1
if "%MODE%"=="" set MODE=debug
if /I "%MODE%"=="debug" (
  set MODEFLAGS=-O0 -g -fcheck=all
) else if /I "%MODE%"=="release" (
  set MODEFLAGS=-O2
) else (
  echo Usage: run_checks.bat [debug^|release]
  exit /b 2
)
set FC=gfortran
set FLAGS=-std=f2018 -ffree-line-length-none -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace %MODEFLAGS%
set BUILD=build\%MODE%
if exist %BUILD% rmdir /s /q %BUILD%
mkdir %BUILD%
set SRC=src\msgarch_kinds.f90 src\msgarch_rng.f90 src\msgarch_special.f90 src\msgarch_distributions.f90 src\msgarch_types.f90 src\msgarch_models.f90 src\msgarch_filter.f90 src\msgarch_simulation.f90 src\msgarch_optimizer.f90 src\msgarch_linalg.f90 src\msgarch_mapping.f90 src\msgarch_estimation.f90 src\msgarch_forecast.f90 src\msgarch_risk.f90 src\msgarch_hmm.f90 src\msgarch_posterior.f90 src\msgarch.f90
%FC% %FLAGS% -J%BUILD% -I%BUILD% %SRC% test\test_helpers.f90 test\test_distributions_models.f90 -o %BUILD%\test_distributions_models.exe || exit /b 1
%FC% %FLAGS% -J%BUILD% -I%BUILD% %SRC% test\test_helpers.f90 test\test_filter_simulation_risk.f90 -o %BUILD%\test_filter_simulation_risk.exe || exit /b 1
%FC% %FLAGS% -J%BUILD% -I%BUILD% %SRC% test\test_helpers.f90 test\test_estimation_hmm.f90 -o %BUILD%\test_estimation_hmm.exe || exit /b 1
%FC% %FLAGS% -J%BUILD% -I%BUILD% %SRC% app\demo_msgarch.f90 -o %BUILD%\demo_msgarch.exe || exit /b 1
%FC% %FLAGS% -J%BUILD% -I%BUILD% %SRC% example\mcmc_example.f90 -o %BUILD%\mcmc_example.exe || exit /b 1
%FC% %FLAGS% -J%BUILD% -I%BUILD% %SRC% example\fit_csv.f90 -o %BUILD%\fit_csv.exe || exit /b 1
%BUILD%\test_distributions_models.exe || exit /b 1
%BUILD%\test_filter_simulation_risk.exe || exit /b 1
%BUILD%\test_estimation_hmm.exe || exit /b 1
%BUILD%\demo_msgarch.exe || exit /b 1
%BUILD%\mcmc_example.exe || exit /b 1
%BUILD%\fit_csv.exe data\returns.csv single sGARCH norm || exit /b 1
%BUILD%\fit_csv.exe data\returns.csv markov sGARCH norm gjrGARCH std || exit /b 1
%BUILD%\fit_csv.exe data\returns.csv mixture sGARCH norm eGARCH ged || exit /b 1
echo %MODE% build, tests, and applications passed.
