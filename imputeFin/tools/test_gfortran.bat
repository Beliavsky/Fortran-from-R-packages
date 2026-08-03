@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build-gfortran
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all
set SRC=imputefin_kinds imputefin_types imputefin_rng imputefin_math imputefin_linalg imputefin_missing imputefin_ar1_gaussian imputefin_ar1_t imputefin_var_t imputefin_wrappers imputefin
for %%S in (%SRC%) do gfortran %FLAGS% -J. -I. -c "%ROOT%\src\%%S.f90" || exit /b 1
for %%T in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% -I. "%%T" *.o -o "%%~nT.exe" || exit /b 1
  "%%~nT.exe" || exit /b 1
)
endlocal
