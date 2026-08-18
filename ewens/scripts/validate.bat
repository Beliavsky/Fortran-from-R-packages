@echo off
setlocal
fpm build || exit /b 1
fpm test || exit /b 1
fpm run --example ewens_demo || exit /b 1
