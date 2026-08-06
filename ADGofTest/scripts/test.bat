@echo off
setlocal
cd /d "%~dp0\.."
if not exist build\check\mod mkdir build\check\mod
if not exist build\check\obj mkdir build\check\obj
if not exist build\check\bin mkdir build\check\bin
set FLAGS=-std=f2018 -pedantic -Wall -Wextra -Werror -O0 -g -fcheck=all -fbacktrace

gfortran %FLAGS% -Jbuild\check\mod -Ibuild\check\mod -c src\adgoftest.f90 -o build\check\obj\adgoftest.o || exit /b 1
for %%F in (test\*.f90) do (
  gfortran %FLAGS% -Jbuild\check\mod -Ibuild\check\mod build\check\obj\adgoftest.o %%F -o build\check\bin\%%~nF.exe || exit /b 1
  build\check\bin\%%~nF.exe || exit /b 1
)
gfortran %FLAGS% -Jbuild\check\mod -Ibuild\check\mod build\check\obj\adgoftest.o example\adgoftest_demo.f90 -o build\check\bin\adgoftest_demo.exe || exit /b 1
build\check\bin\adgoftest_demo.exe || exit /b 1
endlocal
