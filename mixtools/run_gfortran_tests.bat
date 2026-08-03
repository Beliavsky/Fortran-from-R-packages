@echo off
setlocal
set ROOT=%~dp0
set BUILD=%ROOT%build-gfortran
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -O0
gfortran %FLAGS% -c "%ROOT%src\mixtools_kinds.f90" "%ROOT%src\mixtools_status.f90" "%ROOT%src\mixtools_types.f90" "%ROOT%src\mixtools_linalg.f90" "%ROOT%src\mixtools_rng.f90" "%ROOT%src\mixtools_distributions.f90" "%ROOT%src\mixtools_utilities.f90" "%ROOT%src\mixtools_parametric.f90" "%ROOT%src\mixtools_regression.f90" "%ROOT%src\mixtools_semiparametric.f90" "%ROOT%src\mixtools_reliability.f90" "%ROOT%src\mixtools_support.f90" "%ROOT%src\mixtools_diagnostics.f90" "%ROOT%src\mixtools_compat.f90" "%ROOT%src\mixtools.f90"
if errorlevel 1 exit /b 1
for %%T in ("%ROOT%test\*.f90") do (
  gfortran %FLAGS% "%%T" *.o -o "%%~nT.exe"
  if errorlevel 1 exit /b 1
  "%%~nT.exe"
  if errorlevel 1 exit /b 1
)
endlocal
