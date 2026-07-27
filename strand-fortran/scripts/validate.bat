@echo off
setlocal EnableExtensions
cd /d "%~dp0\.."

set MODULES=strand_kinds strand_types strand_linalg strand_stats strand_simplex strand_optimizer strand_data strand_simulation strand
set TESTS=test_data test_exposures_performance test_lifecycle test_optimizer test_simplex test_simulation test_stats
set PROGRAMS=strand_demo portfolio_optimization share_level_simulation

call :build checked "-std=f2018 -O0 -g -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace"
if errorlevel 1 exit /b 1
call :build optimized "-std=f2018 -O2 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror"
if errorlevel 1 exit /b 1
python scripts\audit_release.py
if errorlevel 1 exit /b 1
echo validation: PASS
exit /b 0

:build
set NAME=%~1
set FLAGS=%~2
set DIR=build\validation-%NAME%
if exist "%DIR%" rmdir /s /q "%DIR%"
mkdir "%DIR%\mod" "%DIR%\obj" "%DIR%\bin"
set OBJECTS=
for %%M in (%MODULES%) do (
  gfortran -c %FLAGS% -J "%DIR%\mod" -I "%DIR%\mod" "src\%%M.f90" -o "%DIR%\obj\%%M.o"
  if errorlevel 1 exit /b 1
  call set OBJECTS=%%OBJECTS%% "%DIR%\obj\%%M.o"
)
for %%T in (%TESTS%) do (
  call gfortran %FLAGS% -J "%DIR%\mod" -I "%DIR%\mod" "test\%%T.f90" %%OBJECTS%% -o "%DIR%\bin\%%T.exe"
  if errorlevel 1 exit /b 1
  "%DIR%\bin\%%T.exe"
  if errorlevel 1 exit /b 1
)
for %%P in (strand_demo) do (
  call gfortran %FLAGS% -J "%DIR%\mod" -I "%DIR%\mod" "app\%%P.f90" %%OBJECTS%% -o "%DIR%\bin\%%P.exe"
  if errorlevel 1 exit /b 1
  "%DIR%\bin\%%P.exe"
  if errorlevel 1 exit /b 1
)
for %%P in (portfolio_optimization share_level_simulation) do (
  call gfortran %FLAGS% -J "%DIR%\mod" -I "%DIR%\mod" "example\%%P.f90" %%OBJECTS%% -o "%DIR%\bin\%%P.exe"
  if errorlevel 1 exit /b 1
  "%DIR%\bin\%%P.exe"
  if errorlevel 1 exit /b 1
)
exit /b 0
