@echo off
setlocal
set FC=gfortran
set FLAGS=-std=f2018 -Wall -Wextra -Werror -pedantic -fcheck=all -fbacktrace -O0
if not exist build\checked\mod mkdir build\checked\mod
if not exist build\checked\obj mkdir build\checked\obj
if not exist build\checked\bin mkdir build\checked\bin
%FC% %FLAGS% -Jbuild\checked\mod -Ibuild\checked\mod -c src\ghyp_kinds.f90 src\ghyp_special.f90 src\ghyp_rng.f90 src\ghyp_linalg.f90 src\ghyp_gig.f90 src\ghyp_model.f90 src\ghyp_distribution.f90 src\tsd_types.f90 src\tsd_math.f90 src\tsd_optimize.f90 src\tsd_distributions.f90 src\tsd_fit.f90 src\tsd_moments.f90 src\tsd_spd.f90 src\tsd_profile.f90 src\tsdistributions.f90 src\tsgarch_types.f90 src\tsgarch_model.f90 src\tsgarch_fit.f90 src\tsgarch_simulation.f90 src\tsgarch_forecast.f90 src\tsgarch_diagnostics.f90 src\tsgarch_profile.f90 src\tsgarch_backtest.f90 src\tsgarch_benchmarks.f90 src\tsgarch.f90 || exit /b 1
move /Y *.o build\checked\obj\ >nul
%FC% %FLAGS% -Jbuild\checked\mod -Ibuild\checked\mod -c test\test_support.f90 -o build\checked\obj\test_support.o || exit /b 1
for %%T in (test_models test_simulation_forecast test_estimation test_inference_profile test_backtest test_constraints_vreg) do (
  %FC% %FLAGS% -Jbuild\checked\mod -Ibuild\checked\mod test\%%T.f90 build\checked\obj\*.o -o build\checked\bin\%%T.exe || exit /b 1
  build\checked\bin\%%T.exe || exit /b 1
)
endlocal
