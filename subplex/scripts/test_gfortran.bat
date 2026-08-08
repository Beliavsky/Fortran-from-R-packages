@echo off
setlocal
if "%FC%"=="" set FC=gfortran
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
if exist build-strict rmdir /s /q build-strict
mkdir build-strict
%FC% %FLAGS% -J build-strict -I build-strict -c src\subplex.f90 -o build-strict\subplex.o || exit /b 1
for %%F in (test\*.f90) do (
  %FC% %FLAGS% -J build-strict -I build-strict %%F build-strict\subplex.o -o build-strict\%%~nF.exe || exit /b 1
  build-strict\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90) do (
  %FC% %FLAGS% -J build-strict -I build-strict %%F build-strict\subplex.o -o build-strict\%%~nF.exe || exit /b 1
  build-strict\%%~nF.exe || exit /b 1
)
endlocal
