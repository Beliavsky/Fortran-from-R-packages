@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0
set BUILD=%ROOT%build-gfortran-windows
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -O0 -g -fcheck=all -fbacktrace
set SOURCES=nlme_kinds nlme_status nlme_types nlme_linalg nlme_correlation nlme_variance nlme_pdmat nlme_optimize nlme_covariance nlme_gls nlme_lme nlme_nonlinear nlme_diagnostics nlme_grouped nlme_models nlme segmented_status segmented_types segmented_utils segmented_fit segmented_mixed segmented_inference segmented_wrappers segmented
pushd "%BUILD%\obj"
for %%S in (%SOURCES%) do gfortran %FLAGS% -c "%ROOT%src\%%S.f90" -J "%BUILD%\mod" -I "%BUILD%\mod" || exit /b 1
for %%T in (test_segmented_lm test_stepmented_glm test_inference_selection test_segmented_lme) do (
  gfortran %FLAGS% "%ROOT%test\%%T.f90" *.o -J "%BUILD%\mod" -I "%BUILD%\mod" -o "%BUILD%\bin\%%T.exe" || exit /b 1
  "%BUILD%\bin\%%T.exe" || exit /b 1
)
for %%E in (segmented_linear_example stepmented_example segmented_glm_example segmented_lme_example) do (
  gfortran %FLAGS% "%ROOT%example\%%E.f90" *.o -J "%BUILD%\mod" -I "%BUILD%\mod" -o "%BUILD%\bin\%%E.exe" || exit /b 1
  "%BUILD%\bin\%%E.exe" >nul || exit /b 1
)
gfortran %FLAGS% "%ROOT%app\demo_segmented.f90" *.o -J "%BUILD%\mod" -I "%BUILD%\mod" -o "%BUILD%\bin\demo_segmented.exe" || exit /b 1
"%BUILD%\bin\demo_segmented.exe" >nul || exit /b 1
popd
echo GNU Fortran Windows validation: PASS
