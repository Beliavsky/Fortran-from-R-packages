@echo off
setlocal enabledelayedexpansion
cd /d %~dp0\..
if exist build-strict rmdir /s /q build-strict
mkdir build-strict
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
for %%F in (src\trustoptim_kinds.f90 src\trustoptim_types.f90 src\trustoptim_linalg.f90 src\trustoptim.f90 src\trustoptim_binary.f90) do (
  gfortran %FLAGS% -J build-strict -I build-strict -c %%F -o build-strict\%%~nF.o || exit /b 1
)
set OBJS=build-strict\trustoptim_kinds.o build-strict\trustoptim_types.o build-strict\trustoptim_linalg.o build-strict\trustoptim.o build-strict\trustoptim_binary.o
for %%F in (test\*.f90) do (
  gfortran %FLAGS% -J build-strict -I build-strict %%F %OBJS% -o build-strict\%%~nF.exe || exit /b 1
  build-strict\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90) do (
  gfortran %FLAGS% -J build-strict -I build-strict %%F %OBJS% -o build-strict\%%~nF.exe || exit /b 1
  build-strict\%%~nF.exe || exit /b 1
)
endlocal
