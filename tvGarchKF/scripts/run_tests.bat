@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0\.."
set FC=gfortran
set AR=ar
set FLAGS=-std=f2018 -Wall -Wextra -Werror -pedantic -O0 -g -fcheck=all -fbacktrace
if exist build\windows rmdir /s /q build\windows
mkdir build\windows\mod
mkdir build\windows\obj
mkdir build\windows\bin
set SRC=src\fgarch_kinds.f90 src\fgarch_math.f90 src\fgarch_rng.f90 src\fgarch_distributions.f90 src\fgarch_optimizer.f90 src\fgarch_types.f90 src\fgarch_models.f90 src\fgarch_fit.f90 src\fgarch_risk.f90 src\fgarch.f90 src\tvgarchkf_types.f90 src\tvgarchkf_functions.f90 src\tvgarchkf_core.f90 src\tvgarchkf.f90
%FC% %FLAGS% -J build\windows\mod -I build\windows\mod -c %SRC%
if errorlevel 1 exit /b 1
move /y *.o build\windows\obj\ >nul
%AR% rcs build\windows\libtvgarchkf.a build\windows\obj\*.o
if errorlevel 1 exit /b 1
%FC% %FLAGS% -J build\windows\mod -I build\windows\mod -c test\manual\test_support.f90 -o build\windows\obj\test_support.o
if errorlevel 1 exit /b 1
for %%T in (test_functions test_filter test_simulation test_fit test_tv_parameter) do (
  %FC% %FLAGS% -J build\windows\mod -I build\windows\mod test\manual\%%T.f90 build\windows\obj\test_support.o build\windows\libtvgarchkf.a -o build\windows\bin\%%T.exe
  if errorlevel 1 exit /b 1
  build\windows\bin\%%T.exe
  if errorlevel 1 exit /b 1
)
endlocal
