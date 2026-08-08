@echo off
setlocal
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
if exist build-strict rmdir /s /q build-strict
mkdir build-strict
for %%F in (src\optimflex_types.f90 src\optimflex_linalg.f90 src\optimflex_diff.f90 src\optimflex_helpers.f90 src\optimflex_optimizers.f90 src\optimflex.f90) do (
  gfortran %FLAGS% -c -J build-strict -I build-strict %%F -o build-strict\%%~nF.o || exit /b 1
)
set OBJS=build-strict\optimflex_types.o build-strict\optimflex_linalg.o build-strict\optimflex_diff.o build-strict\optimflex_helpers.o build-strict\optimflex_optimizers.o build-strict\optimflex.o
for %%F in (test\*.f90) do (
  gfortran %FLAGS% -I build-strict %%F %OBJS% -o build-strict\%%~nF.exe || exit /b 1
  build-strict\%%~nF.exe || exit /b 1
)
endlocal
