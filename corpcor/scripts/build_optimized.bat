@echo off
setlocal enabledelayedexpansion
cd /d %~dp0\..
if "%FC%"=="" set FC=gfortran
set BUILD=build\optimized
if exist %BUILD% rmdir /s /q %BUILD%
mkdir %BUILD%\mod %BUILD%\obj %BUILD%\bin
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wno-maybe-uninitialized -O3
set SOURCES=corpcor_kinds corpcor_types corpcor_linalg corpcor_weighted corpcor_matrix_tools corpcor_shrinkage corpcor
set OBJECTS=
for %%S in (%SOURCES%) do (
  %FC% %FLAGS% -J %BUILD%\mod -I %BUILD%\mod -c src\%%S.f90 -o %BUILD%\obj\%%S.o || exit /b 1
  set OBJECTS=!OBJECTS! %BUILD%\obj\%%S.o
)
for %%F in (test\*.f90 example\*.f90 app\*.f90) do (
  %FC% %FLAGS% -I %BUILD%\mod %%F !OBJECTS! -o %BUILD%\bin\%%~nF.exe || exit /b 1
)
for %%E in (%BUILD%\bin\test_*.exe) do %%E || exit /b 1
for %%E in (%BUILD%\bin\example_*.exe) do %%E >nul || exit /b 1
%BUILD%\bin\demo_corpcor.exe >nul || exit /b 1
echo optimized build: PASS
