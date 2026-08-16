@echo off
setlocal
mingw32-make checked
if errorlevel 1 exit /b 1
mingw32-make optimized
