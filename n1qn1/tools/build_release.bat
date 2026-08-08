@echo off
setlocal
if exist build-release rmdir /s /q build-release
mkdir build-release
gfortran -std=f2018 -O3 -Wimplicit-interface -Werror=implicit-interface -J build-release -I build-release src\n1qn1.f90 test\test_n1qn1.f90 -o build-release\test_n1qn1.exe
if errorlevel 1 exit /b 1
build-release\test_n1qn1.exe
if errorlevel 1 exit /b 1
gfortran -std=f2018 -O3 -Wimplicit-interface -Werror=implicit-interface -J build-release -I build-release src\n1qn1.f90 example\banana.f90 -o build-release\banana.exe
if errorlevel 1 exit /b 1
build-release\banana.exe
