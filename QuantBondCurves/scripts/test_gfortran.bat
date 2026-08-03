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
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -ffree-line-length-none -fcheck=all -ffpe-trap=invalid,zero,overflow -O0 -g
set SOURCES=qbc_kinds qbc_status qbc_dates qbc_types qbc_cashflows qbc_curves qbc_bonds qbc_swaps qbc_optimization qbc_calibration quant_bond_curves
set OBJECTS=
for %%U in (%SOURCES%) do (
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
