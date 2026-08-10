@echo off
setlocal EnableExtensions
cd /d "%~dp0\.."
if exist build-strict rmdir /s /q build-strict
mkdir build-strict
set FLAGS=-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all

gfortran %FLAGS% -J build-strict -I build-strict -c src\nnls_kinds.f90 -o build-strict\nnls_kinds.o || exit /b 1
gfortran %FLAGS% -J build-strict -I build-strict -c src\nnls_linalg.f90 -o build-strict\nnls_linalg.o || exit /b 1
gfortran %FLAGS% -J build-strict -I build-strict -c src\nnls_solver.f90 -o build-strict\nnls_solver.o || exit /b 1
gfortran %FLAGS% -J build-strict -I build-strict -c src\nnls.f90 -o build-strict\nnls.o || exit /b 1
set OBJS=build-strict\nnls_kinds.o build-strict\nnls_linalg.o build-strict\nnls_solver.o build-strict\nnls.o

for %%F in (test\*.f90) do (
  gfortran %FLAGS% -I build-strict %OBJS% "%%F" -o "build-strict\%%~nF.exe" || exit /b 1
  "build-strict\%%~nF.exe" || exit /b 1
)
for %%F in (example\*.f90) do (
  gfortran %FLAGS% -I build-strict %OBJS% "%%F" -o "build-strict\%%~nF.exe" || exit /b 1
  "build-strict\%%~nF.exe" || exit /b 1
)
endlocal
