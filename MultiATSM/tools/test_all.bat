@echo off
setlocal EnableExtensions EnableDelayedExpansion

set MODE=%~1
if "%MODE%"=="" set MODE=checked
if /I "%MODE%"=="checked" (
  set FLAGS=-O0 -g -std=f2018 -Wall -Wextra -Werror -fcheck=all -fbacktrace
) else if /I "%MODE%"=="optimized" (
  set FLAGS=-O3 -std=f2018 -Wall -Wextra -Werror
) else (
  echo Usage: %~nx0 [checked^|optimized]
  exit /b 2
)

for %%I in ("%~dp0..") do set ROOT=%%~fI
set BUILD=%ROOT%\build\%MODE%
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%" || exit /b 1
pushd "%BUILD%" || exit /b 1

set SOURCES=multiatsm_kinds multiatsm_linalg multiatsm_types multiatsm_random multiatsm_pca multiatsm_var multiatsm_jll multiatsm_affine multiatsm_likelihood multiatsm_outputs multiatsm_optimization multiatsm_bootstrap multiatsm_bias multiatsm
for %%S in (%SOURCES%) do (
  gfortran %FLAGS% -J. -I. -c "%ROOT%\src\%%S.f90" || goto :fail
)

for %%F in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% -I. "%%F" *.o -llapack -lblas -o "%%~nF.exe" || goto :fail
  "%%~nF.exe" || goto :fail
)

for %%F in ("%ROOT%\example\*.f90" "%ROOT%\app\*.f90") do (
  gfortran %FLAGS% -I. "%%F" *.o -llapack -lblas -o "%%~nF.exe" || goto :fail
  "%%~nF.exe" || goto :fail
)

popd
exit /b 0

:fail
set RC=%ERRORLEVEL%
popd
exit /b %RC%
