@echo off
setlocal enabledelayedexpansion
cd /d %~dp0\..
if exist build-strict rmdir /s /q build-strict
mkdir build-strict
mkdir build-strict\mod
cd build-strict
set FLAGS=-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
for %%F in (bvls_kinds.f90 bvls_types.f90 bvls_qr.f90 bvls_solver.f90 bvls.f90) do (
  gfortran %FLAGS% -Jmod -Imod -c ..\src\%%F || exit /b 1
)
for %%T in (..\test\*.f90) do (
  gfortran %FLAGS% -Imod %%T *.o -o %%~nT.exe || exit /b 1
  %%~nT.exe || exit /b 1
)
for %%E in (..\example\*.f90) do (
  gfortran %FLAGS% -Imod %%E *.o -o %%~nE.exe || exit /b 1
  %%~nE.exe || exit /b 1
)
endlocal
