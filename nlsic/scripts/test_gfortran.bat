@echo off
setlocal EnableExtensions
cd /d "%~dp0\.."
if exist build-strict rmdir /s /q build-strict
mkdir build-strict\mod build-strict\obj build-strict\bin
set F=-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
for %%S in (nlsic_kinds nlsic_types nlsic_linalg nlsic_nnls nlsic_linear nlsic_solver nlsic) do (
  gfortran %F% -J build-strict\mod -I build-strict\mod -c src\%%S.f90 -o build-strict\obj\%%S.o || exit /b 1
)
set O=build-strict\obj\nlsic_kinds.o build-strict\obj\nlsic_types.o build-strict\obj\nlsic_linalg.o build-strict\obj\nlsic_nnls.o build-strict\obj\nlsic_linear.o build-strict\obj\nlsic_solver.o build-strict\obj\nlsic.o
for %%T in (test\*.f90) do (
  gfortran %F% -I build-strict\mod %%T %O% -o build-strict\bin\%%~nT.exe || exit /b 1
  build-strict\bin\%%~nT.exe || exit /b 1
)
for %%T in (example\*.f90) do (
  gfortran %F% -I build-strict\mod %%T %O% -o build-strict\bin\%%~nT.exe || exit /b 1
  build-strict\bin\%%~nT.exe || exit /b 1
)
endlocal
