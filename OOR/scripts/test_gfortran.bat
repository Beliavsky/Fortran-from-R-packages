@echo off
setlocal enabledelayedexpansion
if "%FC%"=="" set FC=gfortran
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
set BUILD=.strict-build
if exist %BUILD% rmdir /s /q %BUILD%
mkdir %BUILD%
for %%F in (src\oor_kinds.f90 src\oor_interfaces.f90 src\oor_random.f90 src\oor_test_functions.f90 src\oor_poo.f90 src\oor_stosoo.f90 src\oor.f90) do (
  %FC% %FLAGS% -c -J %BUILD% -I %BUILD% %%F -o %BUILD%\%%~nF.o || exit /b 1
)
for %%T in (test\*.f90) do (
  %FC% %FLAGS% -J %BUILD% -I %BUILD% %%T %BUILD%\*.o -o %BUILD%\%%~nT.exe || exit /b 1
  %BUILD%\%%~nT.exe || exit /b 1
)
endlocal
