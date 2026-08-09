@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0\.."
if exist build-strict rmdir /s /q build-strict
mkdir build-strict\mod build-strict\obj build-strict\bin
set F=-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
set SRC=nlsr_kinds nlsr_types nlsr_linalg nlsr_derivatives nlsr_core nlsr_models nlsr_stats nlsr
for %%s in (%SRC%) do (
  gfortran %F% -Jbuild-strict\mod -Ibuild-strict\mod -c src\%%s.f90 -o build-strict\obj\%%s.o || exit /b 1
)
set OBJS=
for %%s in (%SRC%) do set OBJS=!OBJS! build-strict\obj\%%s.o
for %%t in (test\*.f90) do (
  gfortran %F% -Ibuild-strict\mod %%t !OBJS! -o build-strict\bin\%%~nt.exe || exit /b 1
  build-strict\bin\%%~nt.exe || exit /b 1
)
for %%e in (example\*.f90) do (
  gfortran %F% -Ibuild-strict\mod %%e !OBJS! -o build-strict\bin\%%~ne.exe || exit /b 1
  build-strict\bin\%%~ne.exe || exit /b 1
)
endlocal
