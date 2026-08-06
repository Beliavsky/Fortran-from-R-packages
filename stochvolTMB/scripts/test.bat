@echo off
setlocal
cd /d %~dp0\..
mingw32-make clean
if errorlevel 1 exit /b 1
mingw32-make MODE=checked test
if errorlevel 1 exit /b 1
mingw32-make clean
if errorlevel 1 exit /b 1
mingw32-make MODE=optimized test
if errorlevel 1 exit /b 1
