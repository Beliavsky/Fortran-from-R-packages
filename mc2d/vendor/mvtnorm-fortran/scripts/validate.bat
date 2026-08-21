@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build-windows
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
set FLAGS=-std=f2018 -pedantic -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace -O0
set SOURCES=mvtnorm_kinds mvtnorm_types mvtnorm_special mvtnorm_linalg mvtnorm_random mvtnorm_distributions mvtnorm_probabilities mvtnorm_triangular mvtnorm_conditioning mvtnorm_quantiles mvtnorm_likelihood mvtnorm
for %%S in (%SOURCES%) do (
  gfortran %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" -c "%ROOT%\src\%%S.f90" -o "%BUILD%\obj\%%S.o" || exit /b 1
)
gfortran %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" -c "%ROOT%\test\test_support.f90" -o "%BUILD%\obj\test_support.o" || exit /b 1
for %%T in (test_conditioning_likelihood test_distributions test_probabilities test_quantiles_scores test_triangular) do (
  gfortran %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" "%ROOT%\test\%%T.f90" "%BUILD%\obj\*.o" -o "%BUILD%\bin\%%T.exe" || exit /b 1
  "%BUILD%\bin\%%T.exe" || exit /b 1
)
echo validation: PASS
