@echo off
setlocal
set FC=gfortran
set FLAGS=-std=f2018 -Wall -Wextra -Werror -fcheck=all -fbacktrace -O0
if exist build\check rmdir /s /q build\check
mkdir build\check
cd build\check
%FC% %FLAGS% -c ..\..\src\rsdc_kinds.f90 ..\..\src\rsdc_types.f90 ..\..\src\rsdc_linalg.f90 ..\..\src\rsdc_random.f90 ..\..\src\rsdc_parameters.f90 ..\..\src\rsdc_filter.f90 ..\..\src\rsdc_likelihood.f90 ..\..\src\vendor\deoptimr_kinds.f90 ..\..\src\vendor\deoptimr_interfaces.f90 ..\..\src\vendor\deoptimr_rng.f90 ..\..\src\vendor\deoptimr_types.f90 ..\..\src\vendor\deoptimr_utils.f90 ..\..\src\vendor\deoptimr_jde.f90 ..\..\src\vendor\deoptimr_ncde.f90 ..\..\src\vendor\deoptimr.f90 ..\..\src\rsdc_simulation.f90 ..\..\src\rsdc_forecast.f90 ..\..\src\rsdc_viterbi.f90 ..\..\src\rsdc_estimation.f90 ..\..\src\rsdc_portfolio.f90 ..\..\src\rsdc_inference.f90 ..\..\src\rsdc_starts.f90 ..\..\src\rsdc_bootstrap.f90 ..\..\src\rsdc_bands.f90 ..\..\src\rsdc.f90
if errorlevel 1 exit /b 1
for %%T in (test_parameters test_filter test_simulation_forecast test_estimation test_portfolio_inference test_statistics) do (
  %FC% %FLAGS% ..\..\test\%%T.f90 *.o -o %%T.exe
  if errorlevel 1 exit /b 1
  %%T.exe
  if errorlevel 1 exit /b 1
)
cd ..\..
endlocal
