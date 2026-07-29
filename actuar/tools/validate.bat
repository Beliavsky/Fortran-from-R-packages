@echo off
setlocal enabledelayedexpansion

set ROOT=%~dp0..
if "%FC%"=="" set FC=gfortran

set DEBUG_FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Wno-compare-reals -fcheck=all -fbacktrace -Werror
set OPT_FLAGS=-std=f2018 -O2 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Wno-compare-reals -Werror

call :run_build debug "%DEBUG_FLAGS%"
if errorlevel 1 exit /b 1
call :run_build optimized "%OPT_FLAGS%"
if errorlevel 1 exit /b 1
exit /b 0

:run_build
set MODE=%~1
set FLAGS=%~2
set BUILD=%ROOT%\build\validation_%MODE%
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
pushd "%BUILD%"

%FC% %FLAGS% -c ^
  "%ROOT%\src\actuar_kinds.f90" ^
  "%ROOT%\src\actuar_special.f90" ^
  "%ROOT%\src\actuar_rng.f90" ^
  "%ROOT%\src\actuar_types.f90" ^
  "%ROOT%\src\actuar_continuous.f90" ^
  "%ROOT%\src\actuar_supplements.f90" ^
  "%ROOT%\src\actuar_discrete.f90" ^
  "%ROOT%\src\actuar_aggregate.f90" ^
  "%ROOT%\src\actuar_phase_type.f90" ^
  "%ROOT%\src\actuar_credibility.f90" ^
  "%ROOT%\src\actuar_grouped.f90" ^
  "%ROOT%\src\actuar_risk.f90" ^
  "%ROOT%\src\actuar.f90"
if errorlevel 1 goto :failed

set OBJECTS=actuar_kinds.o actuar_special.o actuar_rng.o actuar_types.o actuar_continuous.o actuar_supplements.o actuar_discrete.o actuar_aggregate.o actuar_phase_type.o actuar_credibility.o actuar_grouped.o actuar_risk.o actuar.o

for %%T in (test_aggregate test_continuous test_discrete test_grouped_risk test_phase_credibility test_supplements) do (
  %FC% %FLAGS% -I. "%ROOT%\test\%%T.f90" !OBJECTS! -o %%T.exe
  if errorlevel 1 goto :failed
  %%T.exe
  if errorlevel 1 goto :failed
)

%FC% %FLAGS% -I. "%ROOT%\app\actuar_demo.f90" !OBJECTS! -o actuar_demo.exe
if errorlevel 1 goto :failed
actuar_demo.exe >nul
if errorlevel 1 goto :failed

for %%E in (loss_distributions aggregate_and_credibility) do (
  %FC% %FLAGS% -I. "%ROOT%\example\%%E.f90" !OBJECTS! -o %%E.exe
  if errorlevel 1 goto :failed
  %%E.exe >nul
  if errorlevel 1 goto :failed
)

echo validation (%MODE%): PASS
popd
exit /b 0

:failed
popd
exit /b 1
