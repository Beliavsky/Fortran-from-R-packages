@echo off
setlocal
if exist build_gfortran rmdir /s /q build_gfortran
mkdir build_gfortran
cd build_gfortran
set FLAGS=-std=f2018 -Wall -Wextra -Werror -fcheck=all -fbacktrace
set SRC=..\src\fitheavytail_kinds.f90 ..\src\fitheavytail_status.f90 ..\src\fitheavytail_types.f90 ..\src\fitheavytail_linalg.f90 ..\src\fitheavytail_special.f90 ..\src\fitheavytail_rng.f90 ..\src\fitheavytail_tail.f90 ..\src\fitheavytail_elliptical.f90 ..\src\fitheavytail_mvt.f90 ..\src\fitheavytail_mvst.f90 ..\src\fitheavytail.f90 ..\test\test_support.f90
gfortran %FLAGS% -J. -I. -c %SRC% || exit /b 1
for %%T in (test_tail_estimators test_elliptical test_mvt test_mvst) do (
  gfortran %FLAGS% *.o ..\test\%%T.f90 -o %%T.exe || exit /b 1
  %%T.exe || exit /b 1
)
endlocal
