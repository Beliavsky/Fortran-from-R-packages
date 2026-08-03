@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0\.."
if exist build\gfortran-check rmdir /s /q build\gfortran-check
mkdir build\gfortran-check\mod
mkdir build\gfortran-check\bin
set FLAGS=-std=f2018 -Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace -O0
set SOURCES=magic_kinds magic_status magic_tensor magic_square magic_hypercube magic_combinatorics magic
set OBJECTS=
for %%S in (%SOURCES%) do (
  gfortran %FLAGS% -Jbuild\gfortran-check\mod -Ibuild\gfortran-check\mod -c src\%%S.f90 -o build\gfortran-check\%%S.o || exit /b 1
  set OBJECTS=!OBJECTS! build\gfortran-check\%%S.o
)
for %%T in (test\*.f90) do (
  gfortran %FLAGS% -Jbuild\gfortran-check\mod -Ibuild\gfortran-check\mod !OBJECTS! %%T -o build\gfortran-check\bin\%%~nT.exe || exit /b 1
  build\gfortran-check\bin\%%~nT.exe || exit /b 1
)
for %%T in (example\*.f90 app\*.f90) do (
  gfortran %FLAGS% -Jbuild\gfortran-check\mod -Ibuild\gfortran-check\mod !OBJECTS! %%T -o build\gfortran-check\bin\%%~nT.exe || exit /b 1
  build\gfortran-check\bin\%%~nT.exe >nul || exit /b 1
)
endlocal
