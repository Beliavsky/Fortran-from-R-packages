@echo off
mingw32-make clean
if errorlevel 1 exit /b 1
mingw32-make check
if errorlevel 1 exit /b 1
mingw32-make optimized-check
