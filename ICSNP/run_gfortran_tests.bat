@echo off
setlocal
if "%FC%"=="" set FC=gfortran
if "%FFLAGS%"=="" set FFLAGS=-std=f2018 -Wall -Wextra -Werror -pedantic -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -g

if exist build-gfortran rmdir /s /q build-gfortran
mkdir build-gfortran
cd build-gfortran

for %%F in (..\src\icsnp_kinds.f90 ..\src\icsnp_status.f90 ..\src\icsnp_types.f90 ..\src\icsnp_linalg.f90 ..\src\icsnp_special.f90 ..\src\icsnp_pairs.f90 ..\src\icsnp_estimators.f90 ..\src\icsnp_tests.f90 ..\src\icsnp.f90) do (
  %FC% %FFLAGS% -J. -I. -c %%F || exit /b 1
)

set OBJECTS=icsnp_kinds.o icsnp_status.o icsnp_types.o icsnp_linalg.o icsnp_special.o icsnp_pairs.o icsnp_estimators.o icsnp_tests.o icsnp.o

for %%F in (..\test\*.f90) do (
  %FC% %FFLAGS% -J. -I. %OBJECTS% %%F -o %%~nF.exe || exit /b 1
  %%~nF.exe || exit /b 1
)

for %%F in (..\example\*.f90 ..\app\*.f90) do (
  %FC% %FFLAGS% -J. -I. %OBJECTS% %%F -o %%~nF.exe || exit /b 1
  %%~nF.exe || exit /b 1
)
endlocal
