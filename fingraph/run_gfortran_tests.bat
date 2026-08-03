@echo off
setlocal enabledelayedexpansion
set FC=gfortran
set BUILD=build-gfortran
if exist %BUILD% rmdir /s /q %BUILD%
mkdir %BUILD%\mod %BUILD%\obj %BUILD%\bin
set FLAGS=-std=f2008 -pedantic -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all

for %%F in (src\fingraph_kinds.f90 src\fingraph_status.f90 src\fingraph_types.f90 src\fingraph_linalg.f90 src\fingraph_operators.f90 src\fingraph_utils.f90 src\fingraph_learning.f90 src\fingraph_rng.f90 src\fingraph.f90) do (
  %FC% %FLAGS% -J %BUILD%\mod -I %BUILD%\mod -c %%F -o %BUILD%\obj\%%~nF.o || exit /b 1
)
set OBJECTS=
for %%F in (%BUILD%\obj\*.o) do set OBJECTS=!OBJECTS! %%F
for %%F in (test\*.f90) do (
  %FC% %FLAGS% -J %BUILD%\mod -I %BUILD%\mod %%F !OBJECTS! -o %BUILD%\bin\%%~nF.exe || exit /b 1
  %BUILD%\bin\%%~nF.exe || exit /b 1
)
for %%F in (app\*.f90) do (
  %FC% %FLAGS% -J %BUILD%\mod -I %BUILD%\mod %%F !OBJECTS! -o %BUILD%\bin\%%~nF.exe || exit /b 1
  %BUILD%\bin\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90) do (
  %FC% %FLAGS% -J %BUILD%\mod -I %BUILD%\mod %%F !OBJECTS! -o %BUILD%\bin\%%~nF.exe || exit /b 1
  %BUILD%\bin\%%~nF.exe || exit /b 1
)
endlocal
