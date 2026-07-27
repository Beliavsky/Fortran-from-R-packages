@echo off
setlocal
cd /d "%~dp0\.."
if exist build-validation rmdir /s /q build-validation
mkdir build-validation
set FLAGS=-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace -O0 -g -J build-validation -I build-validation
set OBJS=
for %%F in (portfoliooptim_kinds portfoliooptim_types portfoliooptim_linalg portfoliooptim_simplex portfoliooptim_risk portfoliooptim_benders portfoliooptim_projection portfoliooptim) do (
  gfortran %FLAGS% -c src\%%F.f90 -o build-validation\%%F.o || exit /b 1
  call set OBJS=%%OBJS%% build-validation\%%F.o
)
for %%F in (test\test_*.f90) do (
  gfortran %FLAGS% %%F %OBJS% -o build-validation\%%~nF.exe || exit /b 1
  build-validation\%%~nF.exe || exit /b 1
)
echo validation: PASS
endlocal
