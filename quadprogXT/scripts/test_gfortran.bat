@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build-gfortran
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
cd /d "%BUILD%"
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0 -J.

gfortran %FLAGS% -c "%ROOT%\vendor\quadprog-fortran\src\quadprog_kinds.f90" || exit /b 1
gfortran %FLAGS% -c "%ROOT%\vendor\quadprog-fortran\src\quadprog_core.f90" || exit /b 1
gfortran %FLAGS% -c "%ROOT%\vendor\quadprog-fortran\src\quadprog.f90" || exit /b 1
gfortran %FLAGS% -c "%ROOT%\src\quadprogxt.f90" || exit /b 1
gfortran %FLAGS% -c "%ROOT%\test\qpxt_test_support.f90" || exit /b 1

set OBJS=quadprog_kinds.o quadprog_core.o quadprog.o quadprogxt.o qpxt_test_support.o
for %%T in (test_compact_normalize test_base_qp test_absolute_value test_delta_factorized test_full_problem) do (
  gfortran %FLAGS% "%ROOT%\test\%%T.f90" %OBJS% -o %%T.exe || exit /b 1
  %%T.exe || exit /b 1
)

for %%E in (basic_qpxt l1_qpxt) do (
  gfortran %FLAGS% "%ROOT%\example\%%E.f90" quadprog_kinds.o quadprog_core.o quadprog.o quadprogxt.o -o %%E.exe || exit /b 1
  %%E.exe || exit /b 1
)
endlocal
