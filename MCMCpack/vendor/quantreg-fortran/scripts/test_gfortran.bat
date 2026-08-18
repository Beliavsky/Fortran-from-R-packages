@echo off
setlocal
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
if not exist build-strict mkdir build-strict
gfortran %FLAGS% -J build-strict -I build-strict -c src\quantreg_kinds.f90 src\quantreg_types.f90 src\quantreg_linalg.f90 src\quantreg_dense.f90 src\quantreg_select.f90 src\quantreg_local.f90 src\quantreg_nonlinear.f90 src\quantreg_bootstrap.f90 src\quantreg_utils.f90 src\quantreg.f90
if errorlevel 1 exit /b 1
set OBJS=quantreg_kinds.o quantreg_types.o quantreg_linalg.o quantreg_dense.o quantreg_select.o quantreg_local.o quantreg_nonlinear.o quantreg_bootstrap.o quantreg_utils.o quantreg.o
for %%F in (test\*.f90) do (
  gfortran %FLAGS% -I build-strict %%F %OBJS% -o build-strict\%%~nF.exe
  if errorlevel 1 exit /b 1
  build-strict\%%~nF.exe
  if errorlevel 1 exit /b 1
)
for %%F in (example\*.f90) do (
  gfortran %FLAGS% -I build-strict %%F %OBJS% -o build-strict\%%~nF.exe
  if errorlevel 1 exit /b 1
  build-strict\%%~nF.exe
  if errorlevel 1 exit /b 1
)
endlocal
