@echo off
setlocal enabledelayedexpansion
if "%FC%"=="" set FC=gfortran
set ROOT=%~dp0..
set BUILD=%ROOT%\build_gfortran
set MOD=%BUILD%\mod
set OBJ=%BUILD%\obj
set BIN=%BUILD%\bin
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%MOD%" "%OBJ%" "%BIN%"
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -O0 -g
set UNITS=rsolnp_kinds rsolnp_callbacks rsolnp_types rsolnp_linalg rsolnp_evaluate rsolnp_problem rsolnp_solver rsolnp_multistart rsolnp_benchmarks rsolnp
set OBJECTS=
for %%U in (%UNITS%) do (
  %FC% %FLAGS% -J"%MOD%" -I"%MOD%" -c "%ROOT%\src\%%U.f90" -o "%OBJ%\%%U.o" || exit /b 1
  set OBJECTS=!OBJECTS! "%OBJ%\%%U.o"
)
for %%F in ("%ROOT%\test\*.f90") do (
  %FC% %FLAGS% -J"%MOD%" -I"%MOD%" !OBJECTS! "%%F" -o "%BIN%\%%~nF.exe" || exit /b 1
  "%BIN%\%%~nF.exe" || exit /b 1
)
for %%F in ("%ROOT%\example\*.f90" "%ROOT%\app\*.f90") do (
  %FC% %FLAGS% -J"%MOD%" -I"%MOD%" !OBJECTS! "%%F" -o "%BIN%\%%~nF.exe" || exit /b 1
  "%BIN%\%%~nF.exe" || exit /b 1
)
endlocal
