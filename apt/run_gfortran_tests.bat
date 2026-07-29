@echo off
setlocal enabledelayedexpansion
if "%FC%"=="" set FC=gfortran
set COMMON=-std=f2018 -Wall -Wextra -Werror

call :build strict "%COMMON% -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow"
if errorlevel 1 exit /b 1
call :build optimized "%COMMON% -O3"
if errorlevel 1 exit /b 1
echo All apt-fortran builds and tests passed.
exit /b 0

:build
set MODE=%~1
set FLAGS=%~2
set BDIR=build-%MODE%
if exist %BDIR% rmdir /s /q %BDIR%
mkdir %BDIR%\mod %BDIR%\obj %BDIR%\bin
set OBJS=
for %%S in (apt_kinds apt_special apt_regression apt_cointegration apt_ecm apt) do (
  %FC% %FLAGS% -J%BDIR%\mod -I%BDIR%\mod -c src\%%S.f90 -o %BDIR%\obj\%%S.o || exit /b 1
  set OBJS=!OBJS! %BDIR%\obj\%%S.o
)
%FC% %FLAGS% -J%BDIR%\mod -I%BDIR%\mod -c test\test_support.f90 -o %BDIR%\obj\test_support.o || exit /b 1
for %%T in (test_statistics test_reference test_search test_ecm_tests) do (
  %FC% %FLAGS% -J%BDIR%\mod -I%BDIR%\mod !OBJS! %BDIR%\obj\test_support.o test\%%T.f90 -llapack -lblas -o %BDIR%\bin\%%T.exe || exit /b 1
  %BDIR%\bin\%%T.exe || exit /b 1
)
%FC% %FLAGS% -J%BDIR%\mod -I%BDIR%\mod !OBJS! app\apt_demo.f90 -llapack -lblas -o %BDIR%\bin\apt_demo.exe || exit /b 1
%BDIR%\bin\apt_demo.exe > nul || exit /b 1
for %%E in (threshold_search asymmetric_ecm) do (
  %FC% %FLAGS% -J%BDIR%\mod -I%BDIR%\mod !OBJS! example\%%E.f90 -llapack -lblas -o %BDIR%\bin\%%E.exe || exit /b 1
  %BDIR%\bin\%%E.exe > nul || exit /b 1
)
exit /b 0
