@echo off
setlocal
python scripts\audit.py || exit /b 1
fpm build || exit /b 1
fpm test || exit /b 1
fpm run --all || exit /b 1
fpm run --example --all || exit /b 1
