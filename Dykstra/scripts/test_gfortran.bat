@echo off
setlocal
cd /d "%~dp0\.."
if exist build_strict rmdir /s /q build_strict
mkdir build_strict
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0 -J build_strict -I build_strict

gfortran %FLAGS% -c src\dykstra_kinds.f90 -o build_strict\dykstra_kinds.o || exit /b 1
gfortran %FLAGS% -c src\dykstra_linalg.f90 -o build_strict\dykstra_linalg.o || exit /b 1
gfortran %FLAGS% -c src\dykstra_solver.f90 -o build_strict\dykstra_solver.o || exit /b 1
gfortran %FLAGS% -c src\dykstra.f90 -o build_strict\dykstra.o || exit /b 1
set OBJS=build_strict\dykstra_kinds.o build_strict\dykstra_linalg.o build_strict\dykstra_solver.o build_strict\dykstra.o

for %%F in (test\*.f90) do (
    gfortran %FLAGS% "%%F" %OBJS% -o "build_strict\%%~nF.exe" || exit /b 1
    "build_strict\%%~nF.exe" || exit /b 1
)
for %%F in (example\*.f90) do (
    gfortran %FLAGS% "%%F" %OBJS% -o "build_strict\%%~nF.exe" || exit /b 1
    "build_strict\%%~nF.exe" || exit /b 1
)
endlocal
