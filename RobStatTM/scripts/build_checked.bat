@echo off
setlocal
set ROOT=%~dp0..
set BUILD=%ROOT%\build\checked-windows
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
echo Use FPM when available: fpm test --profile debug
fpm test --profile debug
if errorlevel 1 exit /b 1
endlocal
