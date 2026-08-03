@echo off
setlocal
cd /d %~dp0\..
fpm test --profile release --flag "-std=f2018 -O3 -Wall -Wextra -Werror"
if errorlevel 1 exit /b 1
fpm run --example --all --profile release --flag "-std=f2018 -O3 -Wall -Wextra -Werror"
if errorlevel 1 exit /b 1
fpm run --profile release --flag "-std=f2018 -O3 -Wall -Wextra -Werror"
