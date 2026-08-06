@echo off
setlocal
cd /d %~dp0\..
mingw32-make clean
if errorlevel 1 exit /b 1
mingw32-make MODE=checked test example
if errorlevel 1 exit /b 1
mingw32-make MODE=optimized test example
if errorlevel 1 exit /b 1
endlocal
