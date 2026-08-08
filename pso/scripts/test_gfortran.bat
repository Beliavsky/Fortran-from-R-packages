@echo off
setlocal
cd /d "%~dp0\.."
if exist build-gfortran rmdir /s /q build-gfortran
mkdir build-gfortran\mod
mkdir build-gfortran\obj
mkdir build-gfortran\bin
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -Jbuild-gfortran\mod -Ibuild-gfortran\mod

for %%F in (pso_kinds pso_types pso_random pso_lbfgsb pso_core pso_benchmarks pso) do (
  gfortran %FLAGS% -c src\%%F.f90 -o build-gfortran\obj\%%F.o
  if errorlevel 1 exit /b 1
)

set OBJS=build-gfortran\obj\pso_kinds.o build-gfortran\obj\pso_types.o build-gfortran\obj\pso_random.o build-gfortran\obj\pso_lbfgsb.o build-gfortran\obj\pso_core.o build-gfortran\obj\pso_benchmarks.o build-gfortran\obj\pso.o

for %%F in (test\*.f90) do (
  gfortran %FLAGS% %%F %OBJS% -o build-gfortran\bin\%%~nF.exe
  if errorlevel 1 exit /b 1
  build-gfortran\bin\%%~nF.exe
  if errorlevel 1 exit /b 1
)

for %%F in (example\*.f90) do (
  gfortran %FLAGS% %%F %OBJS% -o build-gfortran\bin\%%~nF.exe
  if errorlevel 1 exit /b 1
  build-gfortran\bin\%%~nF.exe
  if errorlevel 1 exit /b 1
)

echo All tests and examples passed.
