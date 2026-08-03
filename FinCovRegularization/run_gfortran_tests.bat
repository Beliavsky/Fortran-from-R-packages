@echo off
setlocal enabledelayedexpansion

if "%FC%"=="" set FC=gfortran
if "%FFLAGS%"=="" set FFLAGS=-std=f2018 -Wall -Wextra -Werror -pedantic -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -O0 -g
set BUILD=build-gfortran

if exist %BUILD% rmdir /s /q %BUILD%
mkdir %BUILD%

set SOURCES=src\fincov_kinds.f90 src\fincov_status.f90 src\fincov_types.f90 src\fincov_utils.f90 src\fincov_rng.f90 src\fincov_linalg.f90 src\fincov_regularization.f90 src\fincov_norms.f90 src\fincov_factor_models.f90 src\fincov_portfolio.f90 src\fincov_cv.f90 src\fincovregularization.f90
set OBJECTS=

for %%S in (%SOURCES%) do (
  %FC% %FFLAGS% -J %BUILD% -I %BUILD% -c %%S -o %BUILD%\%%~nS.o
  if errorlevel 1 exit /b 1
  set OBJECTS=!OBJECTS! %BUILD%\%%~nS.o
)

for %%S in (test\*.f90) do (
  %FC% %FFLAGS% -J %BUILD% -I %BUILD% %%S !OBJECTS! -o %BUILD%\%%~nS.exe
  if errorlevel 1 exit /b 1
  %BUILD%\%%~nS.exe
  if errorlevel 1 exit /b 1
)

for %%S in (app\*.f90 example\*.f90) do (
  %FC% %FFLAGS% -J %BUILD% -I %BUILD% %%S !OBJECTS! -o %BUILD%\%%~nS.exe
  if errorlevel 1 exit /b 1
  %BUILD%\%%~nS.exe
  if errorlevel 1 exit /b 1
)

endlocal
