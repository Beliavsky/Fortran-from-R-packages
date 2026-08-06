@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0\.."

call :build checked "-std=f2018 -O0 -Wall -Wextra -Wconversion-extra -Werror -fcheck=all -fbacktrace"
if errorlevel 1 exit /b 1
call :run checked
if errorlevel 1 exit /b 1

call :build optimized "-std=f2018 -O2 -Wall -Wextra -Wconversion-extra -Werror"
if errorlevel 1 exit /b 1
call :run optimized
if errorlevel 1 exit /b 1

echo validation: PASS
exit /b 0

:build
set MODE=%~1
set FLAGS=%~2
set BUILD=build-%MODE%
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"

for %%S in (ghyp_kinds ghyp_special ghyp_rng ghyp_linalg ghyp_gig ghyp_model ghyp_distribution ghyp_risk ghyp_optimize ghyp_fitting ghyp_portfolio ghyp_utilities ghyp sharpe_rratio_calibration sharpe_rratio_records sharpe_rratio_statistics sharpe_rratio_estimator sharpe_rratio) do (
  gfortran %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" -c "src\%%S.f90" -o "%BUILD%\obj\%%S.o"
  if errorlevel 1 exit /b 1
)

for %%S in (test_records test_calibration test_statistics test_estimator_fixed_nu test_estimator_automatic) do (
  gfortran %FLAGS% -I "%BUILD%\mod" "test\%%S.f90" "%BUILD%\obj\*.o" -o "%BUILD%\bin\%%S.exe"
  if errorlevel 1 exit /b 1
)
gfortran %FLAGS% -I "%BUILD%\mod" "app\sharpe_rratio_demo.f90" "%BUILD%\obj\*.o" -o "%BUILD%\bin\sharpe_rratio_demo.exe"
if errorlevel 1 exit /b 1
gfortran %FLAGS% -I "%BUILD%\mod" "example\known_tail_exponent.f90" "%BUILD%\obj\*.o" -o "%BUILD%\bin\known_tail_exponent.exe"
if errorlevel 1 exit /b 1
exit /b 0

:run
set BUILD=build-%~1
for %%T in (test_records test_calibration test_statistics test_estimator_fixed_nu test_estimator_automatic) do (
  "%BUILD%\bin\%%T.exe"
  if errorlevel 1 exit /b 1
)
"%BUILD%\bin\sharpe_rratio_demo.exe" > nul
if errorlevel 1 exit /b 1
"%BUILD%\bin\known_tail_exponent.exe" > nul
if errorlevel 1 exit /b 1
exit /b 0
