@echo off
setlocal enabledelayedexpansion
set FC=gfortran
set ROOT=%~dp0..
set BUILD=%ROOT%\build_gfortran
set MOD=%BUILD%\mod
set OBJ=%BUILD%\obj
set BIN=%BUILD%\bin

if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%MOD%" "%OBJ%" "%BIN%"

set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -O0 -g
set SOURCES=garchito_kinds garchito_callbacks garchito_types garchito_utils garchito_optimizer garchito_models garchito
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
