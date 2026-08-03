@echo off
setlocal
cd /d "%~dp0.."
fpm clean --all
fpm test --flag "-std=f2018 -Wall -Wextra -Wimplicit-interface -O3"
if errorlevel 1 exit /b 1
fpm run --all --flag "-std=f2018 -Wall -Wextra -Wimplicit-interface -O3"
exit /b %errorlevel%
