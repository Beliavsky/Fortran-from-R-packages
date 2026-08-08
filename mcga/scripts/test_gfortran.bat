@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0\.."
if exist build_strict rmdir /s /q build_strict
mkdir build_strict
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -Jbuild_strict -Ibuild_strict
set OBJS=
for %%F in (src\mcga_kinds.f90 src\mcga_random.f90 src\mcga_bytes.f90 src\mcga_operators.f90 src\mcga_engine.f90 src\mcga.f90) do (
  gfortran %FLAGS% -c %%F -o build_strict\%%~nF.o || exit /b 1
  set OBJS=!OBJS! build_strict\%%~nF.o
)
for %%F in (test\*.f90) do (
  gfortran %FLAGS% %%F !OBJS! -o build_strict\%%~nF.exe || exit /b 1
  build_strict\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90) do (
  gfortran %FLAGS% %%F !OBJS! -o build_strict\%%~nF.exe || exit /b 1
  build_strict\%%~nF.exe || exit /b 1
)
echo All strict tests/examples passed.
