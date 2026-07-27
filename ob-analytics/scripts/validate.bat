@echo off
setlocal
cd /d %~dp0\..
if exist build-validation rmdir /s /q build-validation
mkdir build-validation
cd build-validation
set FLAGS=-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Wno-compare-reals -Werror -fcheck=all -fbacktrace -O0
gfortran %FLAGS% -c -J . -I . ..\src\ob_kinds.f90 ..\src\ob_types.f90 ..\src\ob_utils.f90 ..\src\ob_alignment.f90 ..\src\ob_io.f90 ..\src\ob_events.f90 ..\src\ob_trades.f90 ..\src\ob_depth.f90 ..\src\ob_book.f90 ..\src\ob_processing.f90 ..\src\ob_analytics.f90
if errorlevel 1 exit /b 1
set OBJS=ob_kinds.o ob_types.o ob_utils.o ob_alignment.o ob_io.o ob_events.o ob_trades.o ob_depth.o ob_book.o ob_processing.o ob_analytics.o
for %%F in (..\test\test_*.f90) do (
  gfortran %FLAGS% -I . %%F %OBJS% -o %%~nF.exe
  if errorlevel 1 exit /b 1
)
for %%F in (..\app\*.f90 ..\example\*.f90) do (
  gfortran %FLAGS% -I . %%F %OBJS% -o %%~nF.exe
  if errorlevel 1 exit /b 1
)
cd ..
for %%F in (build-validation\test_*.exe) do (
  %%F
  if errorlevel 1 exit /b 1
)
build-validation\ob_analytics_demo.exe >nul
if errorlevel 1 exit /b 1
build-validation\event_matching.exe >nul
if errorlevel 1 exit /b 1
build-validation\depth_and_book.exe >nul
if errorlevel 1 exit /b 1
echo validation: PASS
