@echo off
setlocal
set FFLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
set BUILD=build_strict
if exist %BUILD% rmdir /s /q %BUILD%
mkdir %BUILD%\mod
mkdir %BUILD%\obj
mkdir %BUILD%\bin
for %%F in (onls_kinds onls_linalg onls_minimize onls_lm onls_core onls) do (
  gfortran %FFLAGS% -J %BUILD%\mod -I %BUILD%\mod -c src\%%F.f90 -o %BUILD%\obj\%%F.o || exit /b 1
)
set OBJS=%BUILD%\obj\onls_kinds.o %BUILD%\obj\onls_linalg.o %BUILD%\obj\onls_minimize.o %BUILD%\obj\onls_lm.o %BUILD%\obj\onls_core.o %BUILD%\obj\onls.o
for %%F in (test\*.f90) do (
  gfortran %FFLAGS% -I %BUILD%\mod %%F %OBJS% -o %BUILD%\bin\%%~nF.exe || exit /b 1
  %BUILD%\bin\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90) do (
  gfortran %FFLAGS% -I %BUILD%\mod %%F %OBJS% -o %BUILD%\bin\%%~nF.exe || exit /b 1
  %BUILD%\bin\%%~nF.exe || exit /b 1
)
echo All strict tests and examples passed.
