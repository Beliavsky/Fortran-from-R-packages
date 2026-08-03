@echo off
setlocal
cd /d "%~dp0.."
fpm clean --all
fpm test --flag "-std=f2018 -Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace -O0 -g"
if errorlevel 1 exit /b 1
fpm run --all --flag "-std=f2018 -Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace -O0 -g"
exit /b %errorlevel%
