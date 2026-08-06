@echo off
setlocal
cd /d %~dp0\..
where fpm >nul 2>nul
if %errorlevel%==0 (
  fpm test --flag "-std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Wconversion-extra -fcheck=all -fbacktrace -ffree-line-length-none"
) else (
  echo FPM was not found. Install Fortran Package Manager and run this script again.
  exit /b 1
)
