@echo off
setlocal enabledelayedexpansion
set MODE=%1
if "%MODE%"=="" set MODE=checked
set ROOT=%~dp0..
set BUILD=%ROOT%\build-%MODE%
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
if /I "%MODE%"=="optimized" (
  set FLAGS=-std=f2018 -O3 -Wall -Wextra -Wpedantic -Wimplicit-interface
) else (
  set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Wimplicit-interface -fcheck=all -fbacktrace
)
for %%S in (pa_kinds pa_types pa_linalg pa_statistics pa_constraints pa_objectives pa_optimizers pa_portfolios pa_views pa_factor_models pa_robust portfolio_analytics) do (
  gfortran !FLAGS! -J "%BUILD%\mod" -I "%BUILD%\mod" -c "%ROOT%\src\%%S.f90" -o "%BUILD%\obj\%%S.o" || exit /b 1
)
gfortran !FLAGS! -J "%BUILD%\mod" -I "%BUILD%\mod" -c "%ROOT%\test\test_support.f90" -o "%BUILD%\obj\test_support.o" || exit /b 1
for %%T in (test_constraints_random test_factor_models test_frontier_rebalancing test_optimization test_statistics test_tail_and_risk_budget test_views) do (
  gfortran !FLAGS! -J "%BUILD%\mod" -I "%BUILD%\mod" "%BUILD%\obj\*.o" "%ROOT%\test\%%T.f90" -o "%BUILD%\bin\%%T.exe" || exit /b 1
  "%BUILD%\bin\%%T.exe" || exit /b 1
)
for %%E in (example_mean_variance example_tail_risk example_black_litterman example_entropy_ranking example_rebalancing) do (
  gfortran !FLAGS! -J "%BUILD%\mod" -I "%BUILD%\mod" "%BUILD%\obj\pa_*.o" "%BUILD%\obj\portfolio_analytics.o" "%ROOT%\example\%%E.f90" -o "%BUILD%\bin\%%E.exe" || exit /b 1
  "%BUILD%\bin\%%E.exe" || exit /b 1
)
gfortran !FLAGS! -J "%BUILD%\mod" -I "%BUILD%\mod" "%BUILD%\obj\pa_*.o" "%BUILD%\obj\portfolio_analytics.o" "%ROOT%\app\demo_portfolio_analytics.f90" -o "%BUILD%\bin\demo_portfolio_analytics.exe" || exit /b 1
"%BUILD%\bin\demo_portfolio_analytics.exe" || exit /b 1
endlocal
