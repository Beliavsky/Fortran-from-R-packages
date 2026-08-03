@echo off
setlocal
call "%~dp0build_backend.bat"
if errorlevel 1 exit /b 1
cd /d "%~dp0.."
if "%~1"=="" (
  fpm run
) else (
  fpm %*
)
exit /b %errorlevel%
