@echo off
setlocal
cd /d %~dp0\..
set SRC=src\optimx_kinds.f90 src\optimx_types.f90 src\optimx_linalg.f90 src\optimx_eval.f90 src\optimx_solvers.f90 src\optimx_checks.f90 src\optimx_api.f90 src\optimx.f90 src\optimx_example_functions.f90
if not exist build mkdir build
gfortran -std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -J build -I build -c %SRC% || exit /b 1
for %%F in (test\*.f90) do (
  gfortran -std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -J build -I build %%F *.o -o build\%%~nF.exe || exit /b 1
  build\%%~nF.exe || exit /b 1
)
del /q *.o 2>nul
endlocal
