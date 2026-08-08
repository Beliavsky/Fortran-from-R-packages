@echo off
setlocal
if "%FC%"=="" set FC=gfortran
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
if exist build_strict rmdir /s /q build_strict
mkdir build_strict
for %%F in (trust_kinds trust_types trust_linalg trust_core trust) do (
  %FC% %FLAGS% -J build_strict -I build_strict -c src\%%F.f90 -o build_strict\%%F.o || exit /b 1
)
set OBJS=build_strict\trust_kinds.o build_strict\trust_types.o build_strict\trust_linalg.o build_strict\trust_core.o build_strict\trust.o
for %%F in (test\*.f90) do (
  %FC% %FLAGS% -J build_strict -I build_strict %%F %OBJS% -o build_strict\%%~nF.exe || exit /b 1
  build_strict\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90) do (
  %FC% %FLAGS% -J build_strict -I build_strict %%F %OBJS% -o build_strict\%%~nF-example.exe || exit /b 1
  build_strict\%%~nF-example.exe || exit /b 1
)
endlocal
