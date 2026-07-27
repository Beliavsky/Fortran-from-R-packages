@echo off
setlocal
cd /d "%~dp0\.."
if exist build rmdir /s /q build
mkdir build\mod build\obj build\bin
set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace -ffree-line-length-132

for %%F in (lsmc_kinds lsmc_math lsmc_random lsmc_linear_algebra lsmc_types lsmc_european lsmc_simulation lsmc_pricing lsmc_surface lsmontecarlo) do (
  gfortran %FLAGS% -Jbuild\mod -Ibuild\mod -c src\%%F.f90 -o build\obj\%%F.o
  if errorlevel 1 exit /b 1
)

set OBJECTS=build\obj\lsmc_kinds.o build\obj\lsmc_math.o build\obj\lsmc_random.o build\obj\lsmc_linear_algebra.o build\obj\lsmc_types.o build\obj\lsmc_european.o build\obj\lsmc_simulation.o build\obj\lsmc_pricing.o build\obj\lsmc_surface.o build\obj\lsmontecarlo.o

for %%F in (test\test_european.f90 test\test_simulation.f90 test\test_american.f90 test\test_exotics_surface.f90 app\lsmontecarlo_demo.f90 example\basic_options.f90 example\quanto_surface.f90) do (
  gfortran %FLAGS% -Jbuild\mod -Ibuild\mod %%F %OBJECTS% -o build\bin\%%~nF.exe
  if errorlevel 1 exit /b 1
  build\bin\%%~nF.exe
  if errorlevel 1 exit /b 1
)

endlocal
