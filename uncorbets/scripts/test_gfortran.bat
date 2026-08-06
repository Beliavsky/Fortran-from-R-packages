@echo off
setlocal
cd /d "%~dp0\.."
if exist build rmdir /s /q build
mingw32-make check
if errorlevel 1 exit /b 1
endlocal
