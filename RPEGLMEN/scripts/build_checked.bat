@echo off
where fpm >nul 2>nul
if %errorlevel%==0 (
  fpm test --flag "-std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Wconversion-extra -fcheck=all -fbacktrace -ffree-line-length-none"
) else (
  mingw32-make MODE=checked clean test
)
