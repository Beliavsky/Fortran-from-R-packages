@echo off
setlocal
cd /d %~dp0\..
if exist build rmdir /s /q build
make test-checked
if errorlevel 1 exit /b 1
make example
