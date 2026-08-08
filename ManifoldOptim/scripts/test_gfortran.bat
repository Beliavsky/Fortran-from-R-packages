@echo off
setlocal
if exist build_gfortran rmdir /s /q build_gfortran
mkdir build_gfortran
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -J build_gfortran -I build_gfortran

gfortran %FLAGS% -c src\manifoldoptim_kinds.f90 src\manifoldoptim_types.f90 src\manifoldoptim_linalg.f90 src\manifoldoptim_manifolds.f90 src\manifoldoptim_solvers.f90 src\manifoldoptim.f90
if errorlevel 1 exit /b 1

for %%F in (test\*.f90) do (
  gfortran %FLAGS% *.o %%F -o build_gfortran\%%~nF.exe
  if errorlevel 1 exit /b 1
  build_gfortran\%%~nF.exe
  if errorlevel 1 exit /b 1
)

del /q *.o 2>nul
endlocal
