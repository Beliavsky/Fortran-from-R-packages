@echo off
fpm test || exit /b 1
fpm run --example basic_workflow || exit /b 1
