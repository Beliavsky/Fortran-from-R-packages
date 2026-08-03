@echo off
setlocal enabledelayedexpansion
cd /d %~dp0\..
if exist build-gfortran rmdir /s /q build-gfortran
mkdir build-gfortran\mod build-gfortran\obj build-gfortran\bin
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
set SOURCES=src\glmnet_kinds.f90 src\glmnet_status.f90 src\glmnet_types.f90 src\glmnet_utils.f90 src\glmnet_gaussian.f90 src\glmnet_glm.f90 src\glmnet_multinomial.f90 src\glmnet_cox.f90 src\glmnet_predict.f90 src\glmnet_assess.f90 src\glmnet_data.f90 src\glmnet_cv.f90 src\glmnet_relax.f90 src\glmnet_control.f90 src\glmnet_api.f90 src\glmnet_families.f90 src\glmnet.f90
for %%F in (%SOURCES%) do (
  gfortran %FLAGS% -Jbuild-gfortran\mod -Ibuild-gfortran\mod -c %%F -o build-gfortran\obj\%%~nF.o || exit /b 1
)
ar rcs build-gfortran\libglmnet.a build-gfortran\obj\*.o || exit /b 1
for %%F in (test\*.f90) do (
  gfortran %FLAGS% -Ibuild-gfortran\mod %%F build-gfortran\libglmnet.a -o build-gfortran\bin\%%~nF.exe || exit /b 1
  build-gfortran\bin\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90 app\*.f90) do (
  gfortran %FLAGS% -Ibuild-gfortran\mod %%F build-gfortran\libglmnet.a -o build-gfortran\bin\%%~nF.exe || exit /b 1
  build-gfortran\bin\%%~nF.exe >nul || exit /b 1
)
echo strict GNU Fortran validation: PASS
