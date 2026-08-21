@echo off
setlocal

if exist build-standalone rmdir /s /q build-standalone
mkdir build-standalone
mkdir build-standalone\mod

gfortran -std=f2018 -Wall -Wextra -Wpedantic -fcheck=all -fbacktrace ^
  -J build-standalone\mod ^
  src\minqa_module.f90 test\test_minqa.f90 ^
  -o build-standalone\test_minqa.exe
if errorlevel 1 exit /b 1

build-standalone\test_minqa.exe
if errorlevel 1 exit /b 1

gfortran -std=f2018 -O3 ^
  -J build-standalone\mod ^
  src\minqa_module.f90 example\minqa_example.f90 ^
  -o build-standalone\minqa_example.exe
if errorlevel 1 exit /b 1

build-standalone\minqa_example.exe
