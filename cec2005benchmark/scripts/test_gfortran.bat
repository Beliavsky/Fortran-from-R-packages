@echo off
setlocal
if "%FC%"=="" set FC=gfortran
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
if exist build-strict rmdir /s /q build-strict
mkdir build-strict\mod
mkdir build-strict\bin
%FC% %FLAGS% -c src\cec2005benchmark.f90 -J build-strict\mod -o build-strict\cec2005benchmark.o || exit /b 1
for %%F in (test\*.f90) do (
  %FC% %FLAGS% -I build-strict\mod %%F build-strict\cec2005benchmark.o -o build-strict\bin\%%~nF.exe || exit /b 1
  build-strict\bin\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90) do (
  %FC% %FLAGS% -I build-strict\mod %%F build-strict\cec2005benchmark.o -o build-strict\bin\example_%%~nF.exe || exit /b 1
  build-strict\bin\example_%%~nF.exe || exit /b 1
)
endlocal
