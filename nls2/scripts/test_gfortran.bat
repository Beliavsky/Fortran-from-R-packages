@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0\.."
if exist build_strict rmdir /s /q build_strict
mkdir build_strict\mod
mkdir build_strict\bin
set FLAGS=-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
gfortran %FLAGS% -c -J build_strict\mod -I build_strict\mod ^
 src\nls2_kinds.f90 src\nls2_types.f90 src\nls2_linalg.f90 src\nls2_random.f90 ^
 src\nls2_core.f90 src\nls2_search.f90 src\nls2_stats.f90 src\nls2.f90 || exit /b 1
move /y *.o build_strict\ >nul
for %%F in (test\test_*.f90) do (
  gfortran %FLAGS% -J build_strict\mod -I build_strict\mod build_strict\*.o "%%F" -o "build_strict\bin\%%~nF.exe" || exit /b 1
  "build_strict\bin\%%~nF.exe" || exit /b 1
)
for %%F in (example\*.f90) do (
  gfortran %FLAGS% -J build_strict\mod -I build_strict\mod build_strict\*.o "%%F" -o "build_strict\bin\%%~nF.exe" || exit /b 1
  "build_strict\bin\%%~nF.exe" || exit /b 1
)
echo All strict tests and examples passed.
