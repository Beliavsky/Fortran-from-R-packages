@echo off
setlocal
if exist build-direct rmdir /s /q build-direct
mkdir build-direct
del /q *.o 2>nul
set FLAGS=-std=f2008 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -O0 -g
gfortran %FLAGS% -J build-direct -I build-direct -c ^
  src\survival_kinds.f90 src\survival_types.f90 src\survival_linalg.f90 src\survival_cox.f90 ^
  src\relsurv_kinds.f90 src\relsurv_ratetable.f90 src\relsurv_parsers.f90 ^
  src\mstate_kinds.f90 src\mstate_types.f90 src\mstate_transitions.f90 src\mstate_data.f90 src\mstate_crprep.f90 ^
  src\mstate_msfit.f90 src\mstate_cox.f90 src\mstate_redrank.f90 src\mstate_probtrans.f90 ^
  src\mstate_nonparametric.f90 src\mstate_simulation.f90 src\mstate_relative.f90 src\mstate_utilities.f90 ^
  src\mstate_markov.f90 src\mstate.f90
if errorlevel 1 exit /b 1
for %%F in (test\*.f90) do (
  gfortran %FLAGS% -J build-direct -I build-direct %%F *.o -o build-direct\%%~nF.exe
  if errorlevel 1 exit /b 1
  build-direct\%%~nF.exe
  if errorlevel 1 exit /b 1
)
for %%F in (example\*.f90) do (
  gfortran %FLAGS% -J build-direct -I build-direct %%F *.o -o build-direct\%%~nF.exe
  if errorlevel 1 exit /b 1
  build-direct\%%~nF.exe
  if errorlevel 1 exit /b 1
)
endlocal
