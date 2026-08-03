@echo off
setlocal
for %%I in ("%~dp0..") do set "ROOT=%%~fI"
call "%ROOT%\scripts\build_backend.bat"
if errorlevel 1 exit /b 1
cd /d "%ROOT%"
fpm %*
