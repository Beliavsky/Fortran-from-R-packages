@echo off
setlocal
cd /d %~dp0\..
fpm test --flag "-std=f2018 -Wall -Wextra -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace"
if errorlevel 1 exit /b 1
fpm run --example --all --flag "-std=f2018 -Wall -Wextra -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace"
if errorlevel 1 exit /b 1
fpm run --flag "-std=f2018 -Wall -Wextra -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace"
