@echo off
setlocal
mingw32-make check
if errorlevel 1 exit /b 1
mingw32-make optimized
if errorlevel 1 exit /b 1
mingw32-make BUILD=build/check MODE=check example
