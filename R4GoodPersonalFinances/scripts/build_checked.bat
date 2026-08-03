@echo off
setlocal
cd /d "%~dp0\.."
if exist build_checked rmdir /s /q build_checked
mkdir build_checked\mod build_checked\obj build_checked\bin
set FLAGS=-std=f2018 -Wall -Wextra -Werror -pedantic -ffree-line-length-none -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
set MODULES=r4gpf_kinds r4gpf_status r4gpf_linalg r4gpf_random r4gpf_optimization r4gpf_mortality r4gpf_finance r4gpf_portfolio r4gpf_household r4gpf_simulation r4good_personal_finances
set OBJS=
for %%M in (%MODULES%) do (
  gfortran %FLAGS% -J build_checked\mod -I build_checked\mod -c src\%%M.f90 -o build_checked\obj\%%M.o || exit /b 1
  call set OBJS=%%OBJS%% build_checked\obj\%%M.o
)
for %%F in (test\*.f90) do (
  gfortran %FLAGS% -I build_checked\mod %OBJS% %%F -o build_checked\bin\%%~nF.exe || exit /b 1
  build_checked\bin\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90 app\*.f90) do (
  gfortran %FLAGS% -I build_checked\mod %OBJS% %%F -o build_checked\bin\%%~nF.exe || exit /b 1
  build_checked\bin\%%~nF.exe || exit /b 1
)
endlocal
