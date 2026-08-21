@echo off
setlocal enabledelayedexpansion
cd /d %~dp0\..
if exist build-validation rmdir /s /q build-validation
mkdir build-validation
set FF=-std=f2008 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -J build-validation -I build-validation
set SRC=src\bzinb_kinds.f90 src\bzinb_special.f90 src\bzinb_rng.f90 src\bzinb_linalg.f90 src\bzinb_optimize.f90 src\bzinb_distributions.f90 src\bzinb_em.f90 src\bzinb_fit.f90 src\bzinb.f90
gfortran %FF% -c %SRC% || exit /b 1
move /y *.o build-validation\ >nul
gfortran %FF% -c test\test_support.f90 -o build-validation\test_support.o || exit /b 1
set OBJS=
for %%F in (build-validation\*.o) do if /I not "%%~nxF"=="test_support.o" set OBJS=!OBJS! %%F
for %%T in (test_distributions test_simulation test_poisson_fits test_bnb_bzinb test_em_parity test_pairwise_full test_weighted_pairwise) do (
  gfortran %FF% test\%%T.f90 build-validation\test_support.o !OBJS! -o build-validation\%%T.exe || exit /b 1
  build-validation\%%T.exe || exit /b 1
)
gfortran %FF% example\demo_bzinb.f90 !OBJS! -o build-validation\demo_bzinb.exe || exit /b 1
build-validation\demo_bzinb.exe || exit /b 1
endlocal
