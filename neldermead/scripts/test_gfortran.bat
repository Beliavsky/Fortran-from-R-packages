@echo off
setlocal enabledelayedexpansion
if "%FC%"=="" set FC=gfortran
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
if exist build-strict rmdir /s /q build-strict
mkdir build-strict\mod build-strict\obj build-strict\bin
set OBJS=
for %%S in (src\neldermead_kinds.f90 src\neldermead_types.f90 src\neldermead_simplex.f90 src\neldermead_core.f90 src\neldermead_frontends.f90 src\neldermead.f90) do (
  %FC% %FLAGS% -Jbuild-strict\mod -Ibuild-strict\mod -c %%S -o build-strict\obj\%%~nS.o || exit /b 1
  set OBJS=!OBJS! build-strict\obj\%%~nS.o
)
for %%T in (test\*.f90) do (
  %FC% %FLAGS% -Jbuild-strict\mod -Ibuild-strict\mod %%T !OBJS! -o build-strict\bin\%%~nT.exe || exit /b 1
  build-strict\bin\%%~nT.exe || exit /b 1
)
for %%E in (example\*.f90) do (
  %FC% %FLAGS% -Jbuild-strict\mod -Ibuild-strict\mod %%E !OBJS! -o build-strict\bin\%%~nE.exe || exit /b 1
  build-strict\bin\%%~nE.exe || exit /b 1
)
