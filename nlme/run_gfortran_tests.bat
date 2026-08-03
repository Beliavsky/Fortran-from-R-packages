@echo off
setlocal enabledelayedexpansion
set MODE=%1
if "%MODE%"=="" set MODE=debug
set ROOT=%~dp0
set BUILD=%ROOT%build-gfortran-%MODE%
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fbacktrace
if /I "%MODE%"=="release" (
  set FLAGS=%FLAGS% -O3
) else (
  set FLAGS=%FLAGS% -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow
)
set SOURCES=nlme_kinds nlme_status nlme_types nlme_linalg nlme_correlation nlme_variance nlme_pdmat nlme_optimize nlme_covariance nlme_gls nlme_lme nlme_nonlinear nlme_diagnostics nlme_grouped nlme_models nlme_test_support nlme
pushd "%BUILD%\obj"
for %%S in (%SOURCES%) do (
  gfortran %FLAGS% -c "%ROOT%src\%%S.f90" -J "%BUILD%\mod" -I "%BUILD%\mod" || exit /b 1
)
for %%T in (test_covariance test_gls test_lme test_nonlinear_diagnostics) do (
  gfortran %FLAGS% "%ROOT%test\%%T.f90" *.o -J "%BUILD%\mod" -I "%BUILD%\mod" -o "%BUILD%\bin\%%T.exe" || exit /b 1
  "%BUILD%\bin\%%T.exe" || exit /b 1
)
popd
echo GNU Fortran %MODE% validation: PASS
