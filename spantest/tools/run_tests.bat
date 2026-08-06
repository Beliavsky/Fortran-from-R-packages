@echo off
setlocal enabledelayedexpansion
set MODE=%1
if "%MODE%"=="" set MODE=checked
if "%FC%"=="" set FC=gfortran
if "%MODE%"=="checked" (
  set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Werror -fcheck=all -fbacktrace
) else if "%MODE%"=="optimized" (
  set FLAGS=-std=f2018 -O3 -Wall -Wextra -Wpedantic -Werror
) else (
  echo usage: run_tests.bat [checked^|optimized]
  exit /b 2
)
set ROOT=%~dp0..
set BUILD=%ROOT%\build\%MODE%
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
for %%S in (spantest_kinds spantest_types spantest_probability spantest_random spantest_linalg spantest_classical spantest_gl spantest_as spantest_simulation spantest) do (
  %FC% %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" -c "%ROOT%\src\%%S.f90" -o "%BUILD%\obj\%%S.o" || exit /b 1
)
for %%T in ("%ROOT%\test\*.f90") do (
  %FC% %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" "%%T" "%BUILD%\obj\*.o" -o "%BUILD%\bin\%%~nT.exe" || exit /b 1
  "%BUILD%\bin\%%~nT.exe" || exit /b 1
)
%FC% %FLAGS% -J"%BUILD%\mod" -I"%BUILD%\mod" "%ROOT%\example\spantest_demo.f90" "%BUILD%\obj\*.o" -o "%BUILD%\bin\spantest_demo.exe" || exit /b 1
"%BUILD%\bin\spantest_demo.exe"
