@echo off
setlocal enabledelayedexpansion
cd /d %~dp0\..
if exist build_strict rmdir /s /q build_strict
mkdir build_strict
set FLAGS=-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
set SRC=cgnm_kinds cgnm_types cgnm_utils cgnm_linalg cgnm_kmeans cgnm_core cgnm_postprocess cgnm_extensions cgnm
set OBJ=
for %%s in (%SRC%) do (
  gfortran %FLAGS% -J build_strict -I build_strict -c src\%%s.f90 -o build_strict\%%s.o || exit /b 1
  set OBJ=!OBJ! build_strict\%%s.o
)
for %%t in (test\*.f90) do (
  gfortran %FLAGS% -I build_strict %%t !OBJ! -o build_strict\%%~nt.exe || exit /b 1
  build_strict\%%~nt.exe || exit /b 1
)
for %%e in (example\*.f90) do (
  gfortran %FLAGS% -I build_strict %%e !OBJ! -o build_strict\%%~ne.exe || exit /b 1
  build_strict\%%~ne.exe || exit /b 1
)
endlocal
