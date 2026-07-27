@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build\windows
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"

set FLAGS=-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace

gfortran -c %FLAGS% ^
 "%ROOT%\src\svdnf_kinds.f90" ^
 "%ROOT%\src\svdnf_types.f90" ^
 "%ROOT%\src\svdnf_stats.f90" ^
 "%ROOT%\src\svdnf_models.f90" ^
 "%ROOT%\src\svdnf_grids.f90" ^
 "%ROOT%\src\svdnf_filter.f90" ^
 "%ROOT%\src\svdnf_simulation.f90" ^
 "%ROOT%\src\svdnf_optimization.f90" ^
 "%ROOT%\src\svdnf.f90" ^
 "%ROOT%\test\test_support.f90" ^
 "%ROOT%\test\test_custom_callbacks.f90"
if errorlevel 1 exit /b 1

set OBJECTS=svdnf_kinds.o svdnf_types.o svdnf_stats.o svdnf_models.o svdnf_grids.o svdnf_filter.o svdnf_simulation.o svdnf_optimization.o svdnf.o
set TEST_OBJECTS=%OBJECTS% test_support.o test_custom_callbacks.o

for %%T in (probabilities models filter jumps forecast optimization) do (
  gfortran %FLAGS% "%ROOT%\test\test_%%T.f90" %TEST_OBJECTS% -o test_%%T.exe
  if errorlevel 1 exit /b 1
  test_%%T.exe
  if errorlevel 1 exit /b 1
)

gfortran %FLAGS% "%ROOT%\app\svdnf_demo.f90" %OBJECTS% -o svdnf_demo.exe
if errorlevel 1 exit /b 1
svdnf_demo.exe > nul
if errorlevel 1 exit /b 1

gfortran %FLAGS% "%ROOT%\example\built_in_models.f90" %OBJECTS% -o built_in_models.exe
if errorlevel 1 exit /b 1
built_in_models.exe > nul
if errorlevel 1 exit /b 1

gfortran %FLAGS% "%ROOT%\example\filter_and_forecast.f90" %OBJECTS% -o filter_and_forecast.exe
if errorlevel 1 exit /b 1
filter_and_forecast.exe > nul
if errorlevel 1 exit /b 1

gfortran %FLAGS% "%ROOT%\example\estimate_taylor.f90" %OBJECTS% -o estimate_taylor.exe
if errorlevel 1 exit /b 1
estimate_taylor.exe > nul
if errorlevel 1 exit /b 1

echo validation: PASS
endlocal
