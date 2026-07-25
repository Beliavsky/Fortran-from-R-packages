@echo off
setlocal enabledelayedexpansion
if not exist build mkdir build
if not exist build\mod mkdir build\mod
if not exist build\obj mkdir build\obj
if not exist build\bin mkdir build\bin
set FLAGS=-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace -Jbuild\mod -Ibuild\mod
set SOURCES=src\fcopulae_kinds.f90 src\fcopulae_rng.f90 src\fcopulae_special.f90 src\fcopulae_integration.f90 src\fcopulae_optimizer.f90 src\fcopulae_linalg.f90 src\fcopulae_distributions.f90 src\fcopulae_utils.f90 src\fcopulae_archimedean.f90 src\fcopulae_elliptical.f90 src\fcopulae_extreme_value.f90 src\fcopulae_empirical.f90 src\fcopulae.f90
for %%F in (%SOURCES%) do (
  gfortran %FLAGS% -c %%F -o build\obj\%%~nF.o || exit /b 1
)
set OBJECTS=build\obj\fcopulae_kinds.o build\obj\fcopulae_rng.o build\obj\fcopulae_special.o build\obj\fcopulae_integration.o build\obj\fcopulae_optimizer.o build\obj\fcopulae_linalg.o build\obj\fcopulae_distributions.o build\obj\fcopulae_utils.o build\obj\fcopulae_archimedean.o build\obj\fcopulae_elliptical.o build\obj\fcopulae_extreme_value.o build\obj\fcopulae_empirical.o build\obj\fcopulae.o
gfortran %FLAGS% app\demo_fcopulae.f90 %OBJECTS% -llapack -lblas -o build\bin\demo_fcopulae.exe || exit /b 1
gfortran %FLAGS% app\fit_csv.f90 %OBJECTS% -llapack -lblas -o build\bin\fit_csv.exe || exit /b 1
gfortran %FLAGS% example\dependence_example.f90 %OBJECTS% -llapack -lblas -o build\bin\dependence_example.exe || exit /b 1
echo Windows release build completed. Tests are run by the Unix build harness used for validation.
