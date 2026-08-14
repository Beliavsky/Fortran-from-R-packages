@echo off
setlocal enabledelayedexpansion
if "%FC%"=="" set FC=gfortran
set ROOT=%~dp0..
set BUILD=%ROOT%\build-validation
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"
set FLAGS=-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -J. -I.
for %%F in (anmc_kinds.f90 anmc_types.f90 anmc_utils.f90 anmc_math.f90 anmc_sampling.f90 anmc_active.f90 anmc_mc.f90 anmc_probabilities.f90 anmc_conservative.f90 anmc.f90) do (
  %FC% %FLAGS% -c "%ROOT%\src\%%F" || exit /b 1
)
for %%T in (test_core test_active test_orthant test_conservative) do (
  %FC% %FLAGS% "%ROOT%\test\%%T.f90" *.o -o %%T.exe || exit /b 1
  %%T.exe || exit /b 1
)
%FC% %FLAGS% "%ROOT%\example\equicorrelated_orthant.f90" *.o -o example_equicorrelated_orthant.exe || exit /b 1
example_equicorrelated_orthant.exe || exit /b 1
endlocal
