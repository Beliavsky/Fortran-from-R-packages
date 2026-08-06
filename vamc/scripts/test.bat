@echo off
setlocal
fpm build
if errorlevel 1 exit /b 1
fpm test
if errorlevel 1 exit /b 1
fpm run
if errorlevel 1 exit /b 1
endlocal
