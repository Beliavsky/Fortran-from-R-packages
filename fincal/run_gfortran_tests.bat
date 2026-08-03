@echo off
setlocal
if "%FC%"=="" set FC=gfortran
if "%FFLAGS%"=="" set FFLAGS=-std=f2018 -Wall -Wextra -Werror -pedantic -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -O0 -g
if exist build-gfortran rmdir /s /q build-gfortran
mkdir build-gfortran

%FC% %FFLAGS% -J build-gfortran -I build-gfortran -c src/fincal_kinds.f90 src/fincal_status.f90 src/fincal_types.f90 src/fincal_rates.f90 src/fincal_tvm.f90 src/fincal_ratios.f90 src/fincal_statistics.f90 src/fincal_accounting.f90 src/fincal.f90
if errorlevel 1 exit /b 1
move /y *.o build-gfortran\ >nul

for %%F in (test\*.f90) do (
  %FC% %FFLAGS% -J build-gfortran -I build-gfortran %%F build-gfortran\*.o -o build-gfortran\%%~nF.exe
  if errorlevel 1 exit /b 1
  build-gfortran\%%~nF.exe
  if errorlevel 1 exit /b 1
)

for %%F in (app\*.f90 example\*.f90) do (
  %FC% %FFLAGS% -J build-gfortran -I build-gfortran %%F build-gfortran\*.o -o build-gfortran\%%~nF.exe
  if errorlevel 1 exit /b 1
  build-gfortran\%%~nF.exe
  if errorlevel 1 exit /b 1
)
endlocal
