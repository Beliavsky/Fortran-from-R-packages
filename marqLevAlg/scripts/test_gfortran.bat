@echo off
setlocal
cd /d %~dp0\..
if exist build-strict rmdir /s /q build-strict
mkdir build-strict
set FF=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
set SRC=src\mla_kinds.f90 src\mla_interfaces.f90 src\mla_linalg.f90 src\mla_derivatives.f90 src\mla_lmm.f90 src\marqlevalg.f90 src\marqlevalg_lmm.f90
for %%F in (test\*.f90) do (
  gfortran %FF% -J build-strict -I build-strict %SRC% %%F -o build-strict\%%~nF.exe || exit /b 1
  build-strict\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90) do (
  gfortran %FF% -J build-strict -I build-strict %SRC% %%F -o build-strict\%%~nF.exe || exit /b 1
  build-strict\%%~nF.exe || exit /b 1
)
endlocal
