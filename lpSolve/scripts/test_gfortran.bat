@echo off
setlocal
cd /d %~dp0\..
if exist build-strict rmdir /s /q build-strict
mkdir build-strict
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
set MODS=-Jbuild-strict -Ibuild-strict

gfortran %FLAGS% %MODS% -c src\lpsolve_types.f90 -o build-strict\lpsolve_types.o || exit /b 1
gfortran %FLAGS% %MODS% -c src\lpsolve_simplex.f90 -o build-strict\lpsolve_simplex.o || exit /b 1
gfortran %FLAGS% %MODS% -c src\lpsolve_core.f90 -o build-strict\lpsolve_core.o || exit /b 1
gfortran %FLAGS% %MODS% -c src\lpsolve_special.f90 -o build-strict\lpsolve_special.o || exit /b 1
gfortran %FLAGS% %MODS% -c src\lpsolve.f90 -o build-strict\lpsolve.o || exit /b 1
set OBJS=build-strict\lpsolve_types.o build-strict\lpsolve_simplex.o build-strict\lpsolve_core.o build-strict\lpsolve_special.o build-strict\lpsolve.o

for %%F in (test\*.f90) do (
  gfortran %FLAGS% %MODS% %%F %OBJS% -o build-strict\%%~nF.exe || exit /b 1
  build-strict\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90) do (
  gfortran %FLAGS% %MODS% %%F %OBJS% -o build-strict\%%~nF.exe || exit /b 1
  build-strict\%%~nF.exe || exit /b 1
)
endlocal
